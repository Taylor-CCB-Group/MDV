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
from .executor import ExecutorConfig, LocalSubprocessExecutor, Handle

logger = logging.getLogger(__name__)


class JobManager:
    """Owner side of the jobs framework. Holds two deliberately separate roots:

    - records_root  — durable owner-side state (ADR-0005), INSIDE the project
      (`<project>/jobs/`), so the catalog never mistakes it for a project and provenance
      travels with the project.
    - workspace_root — ephemeral per-job scratch (ADR-0007), OUTSIDE the project; local
      temp by default, `$SCRATCH`/shared-FS on HPC, wherever the worker computes.

    The job_id is the *only* link between a record and its (possibly remote) workspace.

    `executors` names the executors this manager offers, each with its own workspace root and
    concurrency limit (ADR-0016); the first is the default. Without it, the manager offers one
    executor built from `executor`, `workspace_root` and `max_concurrent`."""

    def __init__(
        self,
        project,
        workspace_root=None,
        records_root=None,
        executor=None,
        max_concurrent: int | None = 2,
        on_submit: Callable[[], None] | None = None,
        executors: dict[str, ExecutorConfig] | None = None,
    ):
        self.project = project
        self.records_root = (
            Path(records_root)
            if records_root is not None
            else Path(project.dir) / JOBS_DIRNAME
        )
        self.executors = executors or {
            "default": ExecutorConfig(
                executor or LocalSubprocessExecutor(max_concurrent),
                Path(workspace_root) if workspace_root is not None else default_workspace_root(project),
                max_concurrent,
            )
        }
        self.default_executor = next(iter(self.executors))
        self.store = JobStore(self.records_root)
        # wake-the-driver hook; the JobService passes its nudge here. None = no-op,
        # so a bare-manager callers (Local/Slurm drivers, tests) work without changing anything
        self.on_submit = on_submit
        self._reconcile_on_boot()  # ADR0005: recover in-flight jobs at startup

    def _executor_name(self, rec: JobRecord) -> str:
        return rec.executor or self.default_executor   # records older than ADR-0016 carry None

    def _config(self, rec: JobRecord) -> ExecutorConfig | None:
        return self.executors.get(self._executor_name(rec))

    def _fail_unconfigured(self, rec: JobRecord) -> None:
        # never hand the record to another executor: its handle means nothing there (ADR-0016)
        self.store.set(
            rec, Status.FAILED,
            error=f"executor {self._executor_name(rec)!r} is not configured on this server",
        )

    def _workspace(self, rec: JobRecord, config: ExecutorConfig) -> Workspace:
        # job_id is the sole correlation key between the durable record (records_root)
        # and the ephemeral scratch (workspace_root) — the two roots are separate by design.
        return Workspace(config.workspace_root, rec.job_id)

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
            config = self._config(rec)
            if config is None:
                self._fail_unconfigured(rec)
            elif rec.handle is not None and config.executor.poll(Handle(**rec.handle)) != "lost":
                self.store.set(rec, Status.RUNNING)                 # reattach
            else:
                self.store.set(rec, Status.QUEUED, handle=None)     # re-queue, since nothing to reattach

    def submit(
        self, tool_id: str, params: dict, executor: str | None = None, rows: list[int] | None = None
    ) -> str:
        spec = get_tool(tool_id)
        validate_params(spec, params, self.project, rows=rows)  # backend-gate
        executor = executor or self.default_executor
        if executor not in self.executors:
            raise ValueError(f"executor {executor!r} is not offered; choose from {list(self.executors)}")
        # write-ahead intent; the executor is stored by name so a later change of default
        # does not move a queued job (ADR-0016)
        rec = self.store.new(tool_id, params, executor=executor, rows=rows)
        if self.on_submit is not None:
            self.on_submit() # wake the driver; it owns dispatch + ingest via tick()
        return rec.job_id


    def _has_slot(self, rec: JobRecord, records: list[JobRecord]) -> bool:
        # the limit is per executor (ADR-0016), so a full local pool holds back no Slurm job.
        # max_concurrent=None -> unbounded: submit every queued job and let the backend's own
        # scheduler queue them (ADR0008); otherwise the manager enforces the bound
        config = self._config(rec)
        if config is None or config.max_concurrent is None:
            return True   # an unconfigured record is let through so dispatch fails it
        name = self._executor_name(rec)
        busy = sum(1 for r in records if r.status in ACTIVE and self._executor_name(r) == name)
        return busy < config.max_concurrent

    def _dispatch(self) -> None:
        # fill every free slot, oldest queued record first among those whose executor has room
        while True:
            records = self.store.load_all()
            nxt = next(
                (r for r in records if r.status == Status.QUEUED.value and self._has_slot(r, records)),
                None,
            )
            if not nxt:
                return
            config = self._config(nxt)
            if config is None:
                self._fail_unconfigured(nxt)
                continue
            spec = get_tool(nxt.tool_id)
            ws = self._workspace(nxt, config)
            self.store.set(nxt, Status.STAGING)
            try:
                self._run_owner_step("materialize", nxt, ws)
                handle = config.executor.submit(
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
            config = self._config(rec)
            if config is None:
                self._fail_unconfigured(rec)
                continue
            ws = self._workspace(rec, config)
            marker = ws.read_marker()
            if marker is None and rec.handle is not None: # marker is primary, if None, check the executor whether the job is still alive
                if config.executor.poll(Handle(**rec.handle)) == "running":
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
