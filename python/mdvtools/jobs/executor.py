from dataclasses import dataclass
from pathlib import Path
from typing import Protocol, Callable
import re
import subprocess
import sys

def _run_cli(argv: list[str]) -> str:
    # Default command-runner: shell out to the Slurm CLI
    return subprocess.run(argv, capture_output=True, text=True, check=True).stdout

def _parse_job_id(sbatch_output: str) -> str:
    # sbatch prints "Submitted job <id>"
    m = re.search(r"Submitted batch job (\d+)", sbatch_output)
    if not m:
        raise ValueError(f"could not parse Slurm job-id from: {sbatch_output!r}")
    return m.group(1)


@dataclass
class Handle:
    kind: str  # 'local' | 'slurm' | 'k8s'
    ref: str  # pid | slurm-id | pod-name


class Executor(Protocol):
    def submit(self, entrypoint: str, workspace: Path) -> Handle: ...
    def poll(self, handle: Handle) -> str: ...  # running | done | lost
    def locate_result(self, handle: Handle, workspace: Path) -> Path: ...


class LocalSubprocessExecutor:
    """ADR0008: local subprocess. Worker dies with owner -> requeue (ADR0005)"""

    def __init__(self, max_concurrent_jobs: int = 2):
        self.max_concurrent_jobs = max_concurrent_jobs
        self._procs: dict[str, subprocess.Popen] = {}

    def submit(self, entrypoint: str, workspace: Path) -> Handle:
        proc = subprocess.Popen(
            [
                sys.executable,
                "-m",
                "mdvtools.jobs.run_worker",
                entrypoint,
                str(workspace),
            ],
            cwd=str(workspace),
        )
        self._procs[str(proc.pid)] = proc
        return Handle("local", str(proc.pid))

    def poll(self, handle: Handle) -> str:
        # Secondary signal; the marker is primary. Catches a worker that vanished without a marker
        proc = self._procs.get(handle.ref)
        if proc is None:
            return "lost"  # no survivor across owner restart -> re-queue
        return "running" if proc.poll() is None else "done"

    def locate_result(self, handle: Handle, workspace: Path) -> Path:
        return workspace / "output"


class SlurmExecutor:
    """
    ADR-0008: slurm executor, submitting jobs using the CLI.
    Precondition (ADR-0010): the compute node shares a POSIX filesystem with the owner at a matching path,
    so the worker reads/writes the same workspace.

    `run` is an injected command-runner (argv -> stdout) for tests mock slurm cli commands
    """
    def __init__(self, run: Callable[[list[str]], str] | None = None, python: str | None = None):
        self._run = run or _run_cli
        self._python = python or sys.executable

    def _render_script(self, entrypoint: str, ws: Path) -> str:
        # Basic directives only for now
        return (
            "#!/bin/bash\n"
            f"#SBATCH --job-name=mdv-{ws.name}\n"
            f"#SBATCH --chdir={ws}\n"
            f"#SBATCH --output={ws / 'slurm-%j.out'}\n"
            f"#SBATCH --error={ws / 'slurm-%j.err'}\n"
            f'{self._python} -m mdvtools.jobs.run_worker "{entrypoint}" "{ws}"'
        )

    def submit(self, entrypoint: str, workspace: Path) -> Handle:
        ws = Path(workspace)
        script_path = ws / "slurm_job.sh"
        script_path.write_text(self._render_script(entrypoint, ws))
        out = self._run(["sbatch", str(script_path)])
        return Handle("slurm", _parse_job_id(out))
