import shutil
from dataclasses import asdict
from pathlib import Path

from .registry import get_tool, validate_params
from .jobstore import JobStore, Status, ACTIVE
from .workspace import Workspace, materialize_columns_tray
from .ingest import ingest_column_output
from .executor import Executor, LocalSubprocessExecutor
from .provenance import build_provenance

# materializers keyed by INPUT shape; ingest handlers keyed by OUTPUT shape
MATERIALIZERS = {"columns": materialize_columns_tray}
INGESTERS = {"column": ingest_column_output}


class JobManager:
    def __init__(self, project, jobs_root, executor=None, max_concurrent=2):
        self.project = project
        self.jobs_root = jobs_root
        self.store = JobStore(self.jobs_root)
        self.executor: Executor = executor or LocalSubprocessExecutor(max_concurrent)
        self.max_concurrent = max_concurrent
        self.store.reconcile_on_boot()  # ADR0005: recover in-flight jobs at startup

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
            ws = Workspace(self.jobs_root, nxt.job_id)
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
            ws = Workspace(self.jobs_root, rec.job_id)
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
