from pathlib import Path
from mdvtools.jobs.executor import Handle, SlurmExecutor

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
