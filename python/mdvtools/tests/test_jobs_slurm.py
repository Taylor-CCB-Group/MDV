from pathlib import Path
from mdvtools.jobs.executor import Handle, SlurmExecutor

def _fake_slurm(squeue="", sacct=""):
    # dispatch a mock Slurm CLI by command name
    def run(argv):
        if argv[0] == "squeue":
            return squeue
        if argv[0] == "sacct":
            return sacct
        raise AssertionError(f"unexpected command {argv!r}")
    return run


def _poll(squeue="", sacct=""):
    return SlurmExecutor(run=_fake_slurm(squeue, sacct)).poll(Handle("slurm", "4242"))


def test_submit_renders_script_and_parses_job_id(tmp_path):
    calls = []

    def fake_run(argv):
        calls.append(argv)
        return "Submitted batch job 4242\n"

    ws = tmp_path / "job-ws"
    ws.mkdir()

    handle = SlurmExecutor(run=fake_run).submit("umap_worker:run", ws)

    # the job-id parsed off sbatch's line, tagged as a slurm handle
    assert handle == Handle("slurm", "4242")

    # sbatch was invoked once, on a generated script file
    assert len(calls) == 1
    assert calls[0][0] == "sbatch"
    script_path = Path(calls[0][1])
    assert script_path.exists()

    body = script_path.read_text()
    assert body.startswith("#!/bin/bash")

    # runs the same sealed run_worker entrypoint against this workspace
    assert "-m mdvtools.jobs.run_worker" in body
    assert "umap_worker:run" in body
    assert str(ws) in body              # --chdir + argv carry the workspace


def test_submit_emits_basic_directives_only_no_resources(tmp_path):
    scripts = []

    def fake_run(argv):
        scripts.append(Path(argv[-1]).read_text())
        return "Submitted batch job 7\n"

    ws = tmp_path / "ws"
    ws.mkdir()
    SlurmExecutor(run=fake_run).submit("umap_worker:run", ws)

    body = scripts[0]
    assert f"#SBATCH --chdir={ws}" in body
    assert "#SBATCH --job-name=" in body
    assert "#SBATCH --output=" in body
    # resource directives are deferred to their own ADR (slice D) — not here yet
    for deferred in ("--mem", "--cpus-per-task", "--time", "--gres"):
        assert deferred not in body


def test_poll_running_while_queued_or_on_node():
    assert _poll(squeue="RUNNING\n") == "running"
    assert _poll(squeue="PENDING\n") == "running"


def test_poll_running_short_circuits_sacct():
    calls = []

    def run(argv):
        calls.append(argv[0])
        return "RUNNING\n" if argv[0] == "squeue" else "COMPLETED\n"

    assert SlurmExecutor(run=run).poll(Handle("slurm", "7")) == "running"
    assert "sacct" not in calls          # squeue was authoritative; no need to hit accounting


def test_poll_done_only_when_accounting_says_completed():
    assert _poll(squeue="", sacct="COMPLETED\n") == "done"


def test_poll_lost_when_scheduler_killed_the_job():
    for state in ("FAILED", "TIMEOUT", "OUT_OF_MEMORY", "CANCELLED", "NODE_FAIL"):
        assert _poll(squeue="", sacct=f"{state}\n") == "lost"


def test_poll_parses_sacct_state_with_suffix():
    assert _poll(squeue="", sacct="CANCELLED by 1000\n") == "lost"


def test_poll_lost_when_no_record_anywhere():
    assert _poll(squeue="", sacct="") == "lost"
