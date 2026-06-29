import shutil
from dataclasses import asdict
from pathlib import Path

from . import JOBS_DIRNAME
from .registry import get_tool, validate_params
from .jobstore import JobStore, Status, ACTIVE
from .workspace import Workspace, materialize_columns_tray, default_workspace_root
from .ingest import ingest_column_output
from .executor import Executor, LocalSubprocessExecutor
from .provenance import build_provenance

# materializers keyed by INPUT shape; ingest handlers keyed by OUTPUT shape
MATERIALIZERS = {"columns": materialize_columns_tray}
INGESTERS = {"column": ingest_column_output}


class JobManager:
    """Owner side of the jobs framework. Holds two deliberately separate roots:

    - records_root  — durable owner-side state (ADR-0005), INSIDE the project
      (`<project>/jobs/`), so the catalog never mistakes it for a project and provenance
      travels with the project.
    - workspace_root — ephemeral per-job scratch (ADR-0007), OUTSIDE the project; local
      temp by default, `$SCRATCH`/shared-FS on HPC, wherever the worker computes.

    The job_id is the *only* link between a record and its (possibly remote) workspace."""

    def __init__(
        self, project, workspace_root=None, records_root=None, executor=None, max_concurrent=2
    ):
        self.project = project
        self.records_root = (
            Path(records_root) if records_root is not None
            else Path(project.dir) / JOBS_DIRNAME
        )
        self.workspace_root = (
            Path(workspace_root) if workspace_root is not None
            else default_workspace_root(project)
        )
        self.store = JobStore(self.records_root)
        self.executor: Executor = executor or LocalSubprocessExecutor(max_concurrent)
        self.max_concurrent = max_concurrent
        self.store.reconcile_on_boot()  # ADR0005: recover in-flight jobs at startup

    def _workspace(self, job_id: str) -> Workspace:
        # job_id is the sole correlation key between the durable record (records_root)
        # and the ephemeral scratch (workspace_root) — the two roots are separate by design.
        return Workspace(self.workspace_root, job_id)

    def submit(self, tool_id: str, params: dict) -> str:
        spec = get_tool(tool_id)
        validate_params(spec, params, self.project)  # backend-gate
        rec = self.store.new(tool_id, params)  # write-ahead intent
        self._dispatch()
        return rec.job_id

    def _busy(self) -> int:
        return sum(1 for r in self.store.load_all() if r.status in ACTIVE)

    def _dispatch(self) -> None:
        # Fill every free slot the executor's bound allows (ADR0008: manager enforces the bound)
        while self._busy() < self.max_concurrent:
            nxt = next(
                (r for r in self.store.load_all() if r.status == Status.QUEUED.value),
                None,
            )
            if not nxt:
                return
            spec = get_tool(nxt.tool_id)
            ws = self._workspace(nxt.job_id)
            self.store.set(nxt, Status.STAGING)
            MATERIALIZERS[spec.input_shape](self.project, spec, nxt.params, ws)
            handle = self.executor.submit(
                spec.entrypoint, ws.root
            )  # submit AFTER intent
            self.store.set(nxt, Status.RUNNING, handle=asdict(handle))

    def tick(self) -> None:
        """Drive periodically. The STATUS marker is the primary completion signal."""
        for rec in self.store.load_all():
            if rec.status != Status.RUNNING.value:
                continue
            ws = self._workspace(rec.job_id)
            marker = ws.read_marker()
            if marker == "done":
                spec = get_tool(rec.tool_id)
                self.store.set(rec, Status.INGESTING)
                manifest = INGESTERS[spec.output.shape](
                    self.project, rec.params, ws
                )  # idempotent
                provenance = build_provenance(rec, manifest)
                self.store.set(rec, Status.DONE, provenance=provenance)
                shutil.rmtree(ws.root, ignore_errors=True)
            elif marker == "failed":
                self.store.set(rec, Status.FAILED)
        self._dispatch()  # a finished job frees up a slot
