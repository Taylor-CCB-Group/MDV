from dataclasses import dataclass
from pathlib import Path
from typing import Protocol
import subprocess
import sys


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
