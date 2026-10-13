from dataclasses import dataclass
from pathlib import Path
from typing import Protocol, Callable
import json
import re
import subprocess
import sys
import urllib.error
import urllib.request

def _run_cli(argv: list[str]) -> str:
    # Default command-runner: shell out to the Slurm CLI
    return subprocess.run(argv, capture_output=True, text=True).stdout

def _parse_job_id(sbatch_output: str) -> str:
    # sbatch prints "Submitted job <id>"
    m = re.search(r"Submitted batch job (\d+)", sbatch_output)
    if not m:
        raise ValueError(f"could not parse Slurm job-id from: {sbatch_output!r}")
    return m.group(1)


def _worker_command(python: str, entrypoint: str, ws: Path) -> str:
    # the line a Slurm job runs, the same for the CLI and the REST executor
    return f'{python} -m mdvtools.jobs.run_worker "{entrypoint}" "{ws}"'


@dataclass
class Handle:
    kind: str  # 'local' | 'slurm' | 'k8s'
    ref: str  # pid | slurm-id | pod-name


class Executor(Protocol):
    def submit(self, entrypoint: str, workspace: Path) -> Handle: ...
    def poll(self, handle: Handle) -> str: ...  # running | done | lost
    def locate_result(self, handle: Handle, workspace: Path) -> Path: ...


@dataclass
class ExecutorConfig:
    """One named executor a manager offers (ADR-0016). Each has its own workspace root, since
    local scratch and the shared Slurm scratch differ; max_concurrent=None = no limit."""
    executor: Executor
    workspace_root: Path
    max_concurrent: int | None


class LocalSubprocessExecutor:
    """ADR0008: local subprocess. Worker dies with owner -> requeue (ADR0005)"""

    def __init__(self, max_concurrent_jobs: int | None = 2):
        self.max_concurrent_jobs = max_concurrent_jobs
        self._procs: dict[str, subprocess.Popen] = {}       # in-memory table of running processes

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
    _ACTIVE_STATES = frozenset({
        "PENDING", "RUNNING", "CONFIGURING", "COMPLETING", "RESIZING", "SUSPENDED", 'REQUEUED'
    })

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
            + _worker_command(self._python, entrypoint, ws)
        )

    def _squeue_state(self, job_id: str) -> str:
        # squeue lists only active jobs; a finished/killed job is absent
        return self._run(["squeue", "-j", job_id, "-h", "-o", "%T"]).strip()

    def _sacct_state(self, job_id: str) -> str:
        # sacct keeps terminal state after the job leaves the queue. -X = allocation only
        # (skip .batch/.extern steps); the state can carry a suffix ("CANCELLED" by 1000)
        out = self._run(["sacct", "-j", job_id, "-n", "-X", "-o", "State"]).strip()
        return out.split()[0] if out else ""

    def submit(self, entrypoint: str, workspace: Path) -> Handle:
        ws = Path(workspace)
        script_path = ws / "slurm_job.sh"
        script_path.write_text(self._render_script(entrypoint, ws))
        out = self._run(["sbatch", str(script_path)])
        return Handle("slurm", _parse_job_id(out))

    def poll(self, handle: Handle) -> str:
        # Slurm only status; the manager will pair this with the STATUS marker in tick()
        if self._squeue_state(handle.ref) in self._ACTIVE_STATES:
            return "running"
        return "done" if self._sacct_state(handle.ref) == "COMPLETED" else "lost"

    def locate_result(self, handle: Handle, workspace: Path) -> Path:
        # shared FS (ADR-0010); the owner reads the same path the worker wrote
        return Path(workspace) / "output"


# (method, path, JSON body or None) -> (HTTP status, JSON body)
RestRequest = Callable[..., tuple[int, dict]]


def rest_request(url: str, token_file: Path | str, user: str | None = None) -> RestRequest:
    """The default slurmrestd transport. The token file is read on every request, so a refresher
    can replace an expiring JWT without restarting the server. Error statuses are returned, not
    raised, so the executor decides what each one means."""

    def request(method: str, path: str, body: dict | None = None) -> tuple[int, dict]:
        headers = {
            "X-SLURM-USER-TOKEN": Path(token_file).read_text().strip(),
            "Content-Type": "application/json",
        }
        if user:
            headers["X-SLURM-USER-NAME"] = user
        data = json.dumps(body).encode() if body is not None else None
        req = urllib.request.Request(url.rstrip("/") + path, data=data, headers=headers, method=method)
        try:
            with urllib.request.urlopen(req) as resp:
                return resp.status, json.loads(resp.read() or b"{}")
        except urllib.error.HTTPError as e:
            try:
                return e.code, json.loads(e.read() or b"{}")
            except ValueError:
                return e.code, {}

    return request


def _slurm_errors(body: dict) -> str:
    return "; ".join(str(e.get("description") or e.get("error")) for e in body.get("errors") or [])


class SlurmRestExecutor:
    """
    ADR-0016: Slurm reached through slurmrestd instead of the CLI, a peer of SlurmExecutor.
    The same shared-filesystem precondition applies (ADR-0010), and handles use the same kind,
    so a job submitted under one mode can be polled under the other.

    `request` is an injected transport ((method, path, body) -> (status, body)); tests fake it.
    """
    _ACTIVE_STATES = SlurmExecutor._ACTIVE_STATES
    # a REST submit does not copy the server's environment; the worker Python is an absolute path
    _JOB_ENVIRONMENT = ["PATH=/usr/bin:/bin"]

    def __init__(self, request: RestRequest, api_version: str = "v0.0.44", python: str | None = None):
        self._request = request
        self._version = api_version
        self._python = python or sys.executable

    def submit(self, entrypoint: str, workspace: Path) -> Handle:
        ws = Path(workspace)
        # slurmrestd does not read #SBATCH lines, so the directives go in the job description
        status, body = self._request("POST", f"/slurm/{self._version}/job/submit", {
            "script": "#!/bin/bash\n" + _worker_command(self._python, entrypoint, ws) + "\n",
            "job": {
                "name": f"mdv-{ws.name}",
                "current_working_directory": str(ws),
                "standard_output": str(ws / "slurm-%j.out"),
                "standard_error": str(ws / "slurm-%j.err"),
                "environment": self._JOB_ENVIRONMENT,
            },
        })
        if status != 200 or "job_id" not in body:
            raise RuntimeError(f"slurmrestd rejected the job ({status}): {_slurm_errors(body)}")
        return Handle("slurm", str(body["job_id"]))

    def poll(self, handle: Handle) -> str:
        # the scheduler holds active jobs, then finished ones for MinJobAge, then 404s
        status, body = self._request("GET", f"/slurm/{self._version}/job/{handle.ref}")
        if status == 200:
            states = {s for job in body.get("jobs") or [] for s in job.get("job_state") or []}
            if states & self._ACTIVE_STATES:
                return "running"
        elif status != 404:
            # for example 511 for an expired token: the job may still be running, so do not lose it
            raise RuntimeError(f"slurmrestd job query failed ({status}): {_slurm_errors(body)}")

        status, body = self._request("GET", f"/slurmdb/{self._version}/job/{handle.ref}")
        if status != 200:
            raise RuntimeError(f"slurmrestd accounting query failed ({status}): {_slurm_errors(body)}")
        states = {s for job in body.get("jobs") or [] for s in (job.get("state") or {}).get("current") or []}
        return "done" if "COMPLETED" in states else "lost"

    def locate_result(self, handle: Handle, workspace: Path) -> Path:
        # shared FS (ADR-0010); the owner reads the same path the worker wrote
        return Path(workspace) / "output"
