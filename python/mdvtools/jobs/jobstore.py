from dataclasses import dataclass, asdict, field
import dataclasses
from enum import Enum
from pathlib import Path
import json
import time
import uuid
import os
import logging

import numpy as np

from .provenance import rows_hash

logger = logging.getLogger(__name__)


class Status(str, Enum):
    QUEUED = "queued"
    STAGING = "staging"
    RUNNING = "running"
    INGESTING = "ingesting"
    DONE = "done"
    FAILED = "failed"
    CANCELLED = "cancelled"
    STALE = "stale"
    LOST = "lost"


@dataclass
class JobRecord:
    job_id: str
    tool_id: str
    params: dict
    status: str = Status.QUEUED.value
    handle: dict | None = (
        None  # None until submitted (the submit <-> record race window)
    )
    input_filter_hash: str | None = (
        None  # subset-taking tools only. concat_columns -> None
    )
    created: float = field(default_factory=time.time)
    provenance: dict | None = None  # promoted at ingest (ADR0007)
    error: str | None = None # set when the owner side fails this record
    executor: str | None = None  # name chosen at submit (ADR-0016); None = the manager's default
    selection_rows: int | None = None  # row count of the selection (ADR-0017); None = all rows


# states that mean "work was in flight when we stopped" - recoverable not terminal
ACTIVE = (Status.STAGING.value, Status.RUNNING.value, Status.INGESTING.value)


class JobStore:
    """Durable owner side (ADR0005). One JSON per job; the in-memory record is a cache.

    Records live *inside* the project (`<project>/jobs/records/`) — separate from the
    ephemeral per-job scratch (ADR-0007), which lives outside the project. The two are
    linked only by job_id."""

    def __init__(self, records_root: Path):
        self.records_dir = Path(records_root) / "records"
        self.records_dir.mkdir(parents=True, exist_ok=True)
        # malformed records moved aside by load_all; outside the *.json glob of records_dir (ADR0012)
        self.quarantine_dir = self.records_dir / "quarantine"
        # a selection's row indices, one sidecar per job, kept as long as its record (ADR-0017)
        self.selections_dir = Path(records_root) / "selections"

    def new(
        self, tool_id: str, params: dict, executor: str | None = None, rows: list[int] | None = None
    ) -> JobRecord:
        rec = JobRecord(
            uuid.uuid4().hex[:12], tool_id, params, executor=executor
        )  # job_id fixed at submit
        if rows is not None:
            # sidecar first, so a record never points at rows that are not on disk
            selected = np.sort(np.asarray(rows, dtype=np.uint32))
            self.selections_dir.mkdir(exist_ok=True)
            path = self.rows_path(rec.job_id)
            tmp = path.with_suffix(".npy.tmp")
            with open(tmp, "wb") as f:
                np.save(f, selected)
            os.replace(tmp, path)
            rec.input_filter_hash = rows_hash(selected)
            rec.selection_rows = len(selected)
        self._write(rec)
        return rec

    def rows_path(self, job_id: str) -> Path:
        return self.selections_dir / f"{job_id}.npy"

    def load_rows(self, job_id: str) -> np.ndarray | None:
        """The sorted row indices a job runs on, or None when it runs on all rows."""
        path = self.rows_path(job_id)
        return np.load(path) if path.exists() else None

    def _write(self, rec: JobRecord) -> None:
        path = self.records_dir / f"{rec.job_id}.json"
        tmp = path.with_suffix(".json.tmp")
        tmp.write_text(json.dumps(asdict(rec)))
        os.replace(tmp, path) # atomic rename: readers see the old record or the new one

    def set(self, rec: JobRecord, status: Status, **fields) -> JobRecord:
        rec.status = status.value
        valid = {f.name for f in dataclasses.fields(rec)}
        for k, v in fields.items():
            if k not in valid:
                raise AttributeError(f"Unknown JobRecord field {k!r}")
            setattr(rec, k, v)
        self._write(rec)
        return rec

    def load_all(self) -> list[JobRecord]:
        """Parse each record on its own (ADR0012)."""
        records = []
        for p in self.records_dir.glob("*.json"):
            try:
                records.append(JobRecord(**json.loads(p.read_text())))
            except (ValueError, TypeError):
                # terminal: bad JSON or wrong fields will not heal, so move it aside for inspection
                self.quarantine_dir.mkdir(exist_ok=True)
                os.replace(p, self.quarantine_dir / p.name)
                logger.exception("job record %s is malformed; quarantined", p.name)
            except OSError:
                # transient: a lock or IO hiccup, so leave it in place for the next pass
                logger.warning("job record %s unreadable this pass; will retry", p.name, exc_info=True)
        return records
