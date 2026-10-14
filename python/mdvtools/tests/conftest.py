import json
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

import pytest


class _FakeSlurmrestdHandler(BaseHTTPRequestHandler):
    """Records each request; a POST submits job 1, a GET finds no job (404)."""
    seen: list = []

    def _reply(self, status, body):
        payload = json.dumps(body).encode()
        self.send_response(status)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(payload)))
        self.end_headers()
        self.wfile.write(payload)

    def do_POST(self):
        length = int(self.headers.get("Content-Length", 0))
        self.seen.append({
            "path": self.path,
            "token": self.headers.get("X-SLURM-USER-TOKEN"),
            "user": self.headers.get("X-SLURM-USER-NAME"),
            "body": json.loads(self.rfile.read(length)),
        })
        self._reply(200, {"job_id": 1, "errors": []})

    def do_GET(self):
        self._reply(404, {"jobs": [], "errors": [{"description": "gone"}]})

    def log_message(self, format, *args):
        pass    # keep request logs out of the test output


@pytest.fixture
def fake_slurmrestd():
    """A local HTTP server standing in for slurmrestd; yields (url, recorded requests)."""
    seen: list = []
    handler = type("Handler", (_FakeSlurmrestdHandler,), {"seen": seen})
    server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    yield f"http://127.0.0.1:{server.server_address[1]}", seen
    server.shutdown()
    server.server_close()
