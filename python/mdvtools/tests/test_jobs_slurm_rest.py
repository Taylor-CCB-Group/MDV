import json
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

import pytest

from mdvtools.jobs.executor import Handle, SlurmRestExecutor, rest_request


def _fake_api(routes):
    """A fake slurmrestd: maps (method, path) to (status, body) and records each call."""
    calls = []

    def request(method, path, body=None):
        calls.append((method, path, body))
        if (method, path) not in routes:
            raise AssertionError(f"unexpected request {method} {path}")
        return routes[(method, path)]

    return request, calls


def _scheduler(*states):
    return 200, {"jobs": [{"job_state": list(states)}], "errors": []}


def _accounting(*states):
    return 200, {"jobs": [{"state": {"current": list(states)}}], "errors": []}


_GONE = (404, {"jobs": [], "errors": [{"error_number": 2017, "description": "Unable to query JobId=7"}]})


def _poll(scheduler, accounting=None):
    routes = {("GET", "/slurm/v0.0.44/job/7"): scheduler}
    if accounting is not None:
        routes[("GET", "/slurmdb/v0.0.44/job/7")] = accounting
    request, _ = _fake_api(routes)
    return SlurmRestExecutor(request).poll(Handle("slurm", "7"))


def test_submit_posts_the_worker_job_and_returns_a_slurm_handle(tmp_path):
    request, calls = _fake_api({
        ("POST", "/slurm/v0.0.42/job/submit"): (200, {"job_id": 4242, "errors": []}),
    })
    ws = tmp_path / "job-ws"
    ws.mkdir()

    handle = SlurmRestExecutor(request, api_version="v0.0.42", python="/shared/venv/bin/python").submit(
        "umap_worker:run", ws
    )

    # the same handle kind as the CLI executor, so jobs survive a switch between cli and rest
    assert handle == Handle("slurm", "4242")
    [(_, _, body)] = calls
    assert body["script"].startswith("#!/bin/bash")
    assert f'/shared/venv/bin/python -m mdvtools.jobs.run_worker "umap_worker:run" "{ws}"' in body["script"]
    job = body["job"]
    assert job["current_working_directory"] == str(ws)
    assert job["standard_output"] == str(ws / "slurm-%j.out")
    assert job["standard_error"] == str(ws / "slurm-%j.err")
    assert job["name"] == "mdv-job-ws"
    # REST does not copy the server's environment, and server secrets must not reach Slurm
    assert job["environment"] == ["PATH=/usr/bin:/bin"]


def test_submit_raises_slurms_message_when_rejected(tmp_path):
    request, _ = _fake_api({
        ("POST", "/slurm/v0.0.44/job/submit"): (
            500, {"errors": [{"error_number": 2008, "description": "Invalid partition name specified"}]},
        ),
    })

    with pytest.raises(RuntimeError, match="Invalid partition name specified"):
        SlurmRestExecutor(request).submit("umap_worker:run", tmp_path)


def test_poll_running_while_the_scheduler_holds_an_active_job():
    assert _poll(_scheduler("RUNNING")) == "running"
    assert _poll(_scheduler("PENDING")) == "running"


def test_poll_running_does_not_ask_accounting():
    request, calls = _fake_api({("GET", "/slurm/v0.0.44/job/7"): _scheduler("RUNNING")})

    SlurmRestExecutor(request).poll(Handle("slurm", "7"))

    assert [path for _, path, _ in calls] == ["/slurm/v0.0.44/job/7"]


def test_poll_done_only_when_accounting_says_completed():
    # the scheduler keeps a finished job for MinJobAge, then forgets it
    assert _poll(_scheduler("COMPLETED"), _accounting("COMPLETED")) == "done"
    assert _poll(_GONE, _accounting("COMPLETED")) == "done"


def test_poll_lost_when_slurm_killed_the_job_or_has_no_record():
    for state in ("FAILED", "TIMEOUT", "OUT_OF_MEMORY", "CANCELLED", "NODE_FAIL"):
        assert _poll(_GONE, _accounting(state)) == "lost"
    assert _poll(_GONE, (200, {"jobs": [], "errors": []})) == "lost"


def test_poll_raises_on_an_unexpected_status_rather_than_losing_the_job():
    # slurmrestd answers 511 to an expired token; the job may still be running
    with pytest.raises(RuntimeError, match="511"):
        _poll((511, {"jobs": [], "errors": [{"error_number": 1007, "description": "Unable to query JobId=7"}]}))


def test_locate_result_is_workspace_output(tmp_path):
    executor = SlurmRestExecutor(lambda method, path, body=None: (200, {}))
    assert executor.locate_result(Handle("slurm", "7"), tmp_path) == tmp_path / "output"


class _RecordingHandler(BaseHTTPRequestHandler):
    seen: list = []

    def do_POST(self):
        length = int(self.headers.get("Content-Length", 0))
        _RecordingHandler.seen.append({
            "path": self.path,
            "token": self.headers.get("X-SLURM-USER-TOKEN"),
            "user": self.headers.get("X-SLURM-USER-NAME"),
            "body": json.loads(self.rfile.read(length)),
        })
        payload = json.dumps({"job_id": 1, "errors": []}).encode()
        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(payload)))
        self.end_headers()
        self.wfile.write(payload)

    def do_GET(self):
        payload = json.dumps({"jobs": [], "errors": [{"description": "gone"}]}).encode()
        self.send_response(404)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(payload)))
        self.end_headers()
        self.wfile.write(payload)

    def log_message(self, format, *args):
        pass    # keep request logs out of the test output


@pytest.fixture
def fake_slurmrestd():
    _RecordingHandler.seen = []
    server = ThreadingHTTPServer(("127.0.0.1", 0), _RecordingHandler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    yield f"http://127.0.0.1:{server.server_address[1]}"
    server.shutdown()
    server.server_close()


def test_rest_request_reads_the_token_file_on_every_request(tmp_path, fake_slurmrestd):
    token_file = tmp_path / "token"
    token_file.write_text("first\n")
    request = rest_request(fake_slurmrestd, token_file, user="mdv")

    assert request("POST", "/slurm/v0.0.44/job/submit", {"script": "x"}) == (200, {"job_id": 1, "errors": []})
    token_file.write_text("second\n")    # a refresher replaced the expiring token
    request("POST", "/slurm/v0.0.44/job/submit", {"script": "y"})

    assert [s["token"] for s in _RecordingHandler.seen] == ["first", "second"]
    assert _RecordingHandler.seen[0]["user"] == "mdv"
    assert _RecordingHandler.seen[0]["path"] == "/slurm/v0.0.44/job/submit"
    assert _RecordingHandler.seen[0]["body"] == {"script": "x"}


def test_rest_request_returns_error_statuses_instead_of_raising(tmp_path, fake_slurmrestd):
    token_file = tmp_path / "token"
    token_file.write_text("t")

    status, body = rest_request(fake_slurmrestd, token_file)("GET", "/slurm/v0.0.44/job/7")

    assert status == 404
    assert body["errors"] == [{"description": "gone"}]
