from dataclasses import dataclass, asdict, field
import dataclasses
from enum import Enum
from pathlib import Path
import json
import time
import uuid
import os
import logging

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

    def new(self, tool_id: str, params: dict) -> JobRecord:
        rec = JobRecord(
            uuid.uuid4().hex[:12], tool_id, params
        )  # job_id fixed at submit
        self._write(rec)
        return rec

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
