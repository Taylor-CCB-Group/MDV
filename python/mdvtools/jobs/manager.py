import json
import shutil
import subprocess
import sys
from dataclasses import asdict
from collections.abc import Callable
from pathlib import Path
import logging

from . import JOBS_DIRNAME
from .registry import get_tool, validate_params
from .jobstore import JobRecord, JobStore, Status, ACTIVE
from .workspace import Workspace, default_workspace_root
from .executor import Executor, LocalSubprocessExecutor, Handle

logger = logging.getLogger(__name__)


class JobManager:
    """Owner side of the jobs framework. Holds two deliberately separate roots:

    - records_root  — durable owner-side state (ADR-0005), INSIDE the project
      (`<project>/jobs/`), so the catalog never mistakes it for a project and provenance
      travels with the project.
    - workspace_root — ephemeral per-job scratch (ADR-0007), OUTSIDE the project; local
      temp by default, `$SCRATCH`/shared-FS on HPC, wherever the worker computes.

    The job_id is the *only* link between a record and its (possibly remote) workspace."""

    def __init__(
        self,
        project,
        workspace_root=None,
        records_root=None,
        executor=None,
        max_concurrent: int | None = 2,
        on_submit: Callable[[], None] | None = None,
    ):
        self.project = project
        self.records_root = (
            Path(records_root)
            if records_root is not None
            else Path(project.dir) / JOBS_DIRNAME
        )
        self.workspace_root = (
            Path(workspace_root)
            if workspace_root is not None
            else default_workspace_root(project)
        )
        self.store = JobStore(self.records_root)
        self.executor: Executor = executor or LocalSubprocessExecutor(max_concurrent)
        self.max_concurrent = max_concurrent
        # wake-the-driver hook; the JobService passes its nudge here. None = no-op,
        # so a bare-manager callers (Local/Slurm drivers, tests) work without changing anything
        self.on_submit = on_submit
        self._reconcile_on_boot()  # ADR0005: recover in-flight jobs at startup

    def _workspace(self, job_id: str) -> Workspace:
        # job_id is the sole correlation key between the durable record (records_root)
        # and the ephemeral scratch (workspace_root) — the two roots are separate by design.
        return Workspace(self.workspace_root, job_id)

    def _run_owner_step(self, step: str, rec: JobRecord, ws: Workspace) -> str:
        """Run an owner step in its own process and wait for it (ADR-0014): the HDF5 work stays off
        the server's gevent hub and out of the process its request handlers use. Returns the
        step's stdout; a failed step raises with the last line of its stderr."""
        result = subprocess.run(
            [sys.executable, "-m", "mdvtools.jobs.owner_step", step, self.project.dir, str(ws.root)],
            input=json.dumps(asdict(rec)),
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            logger.error("job %s: %s step failed:\n%s", rec.job_id, step, result.stderr)
            lines = result.stderr.strip().splitlines()
            raise RuntimeError(lines[-1] if lines else f"{step} step exited with {result.returncode}")
        return result.stdout

    def _reconcile_on_boot(self) -> None:
        # rather than blindly re-queuing every active record, poll its durable handle
        # and only re-queue if the handle is not present or the job is not running
        for rec in self.store.load_all():
            if rec.status not in ACTIVE:
                continue
            if rec.handle is not None and self.executor.poll(Handle(**rec.handle)) != "lost":
                self.store.set(rec, Status.RUNNING)                 # reattach
            else:
                self.store.set(rec, Status.QUEUED, handle=None)     # re-queue, since nothing to reattach

    def submit(self, tool_id: str, params: dict) -> str:
        spec = get_tool(tool_id)
        validate_params(spec, params, self.project)  # backend-gate
        rec = self.store.new(tool_id, params)  # write-ahead intent
        if self.on_submit is not None:
            self.on_submit() # wake the driver; it owns dispatch + ingest via tick()
        return rec.job_id


    def _busy(self) -> int:
        return sum(1 for r in self.store.load_all() if r.status in ACTIVE)

    def _dispatch(self) -> None:
        # max_concurrent=None -> unbounded: submit every queued job and let the backend's own
        # scheduler queue them (ADR0008) otherwise
        # fill every free slot the executor's bound allows (ADR0008: manager enforces the bound for LocalSubprocessExecutor)
        while self.max_concurrent is None or self._busy() < self.max_concurrent:
            nxt = next(
                (r for r in self.store.load_all() if r.status == Status.QUEUED.value),
                None,
            )
            if not nxt:
                return
            spec = get_tool(nxt.tool_id)
            ws = self._workspace(nxt.job_id)
            self.store.set(nxt, Status.STAGING)
            try:
                self._run_owner_step("materialize", nxt, ws)
                handle = self.executor.submit(
                    spec.entrypoint, ws.root
                ) # submit AFTER intent
            except Exception as e:
                # owner side failure: fail this record, free its slot, try the next queued (ADR:0012)
                logger.exception("job %s: dispatch failed", nxt.job_id)
                self.store.set(nxt, Status.FAILED, error=str(e))
                continue
            self.store.set(nxt, Status.RUNNING, handle=asdict(handle))

    def tick(self) -> None:
        """Drive periodically. The STATUS marker is the primary completion signal."""
        for rec in self.store.load_all():
            if rec.status != Status.RUNNING.value:
                continue
            ws = self._workspace(rec.job_id)
            marker = ws.read_marker()
            if marker is None and rec.handle is not None: # marker is primary, if None, check the executor whether the job is still alive
                if self.executor.poll(Handle(**rec.handle)) == "running":
                    continue
                marker = ws.read_marker() or "failed"
            if marker == "done":
                try:
                    self.store.set(rec, Status.INGESTING)
                    # idempotent: writes the columns and their provenance pointers
                    provenance = json.loads(self._run_owner_step("ingest", rec, ws))
                    # commit last: pointer + column data are written while status is still "INGESTING"
                    self.store.set(rec, Status.DONE, provenance=provenance)
                    shutil.rmtree(ws.root, ignore_errors=True)
                except Exception as e:
                    # owner side failure: fail this record and carry on with the rest (ADR:0012)
                    logger.exception("job %s: ingest failed", rec.job_id)
                    self.store.set(rec, Status.FAILED, error=str(e))
            elif marker == "failed":
                self.store.set(rec, Status.FAILED)
        self._dispatch()  # a finished job frees up a slot
