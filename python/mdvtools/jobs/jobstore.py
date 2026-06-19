from dataclasses import dataclass, asdict, field
import dataclasses
from enum import Enum
from pathlib import Path
import json
import time
import uuid


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


# states that mean "work was in flight when we stopped" - recoverable not terminal
ACTIVE = (Status.STAGING.value, Status.RUNNING.value, Status.INGESTING.value)


class JobStore:
    """Durable owner side (ADR0005). One JSON per job; the in-memory record is a cache."""

    def __init__(self, jobs_root: Path):
        self.records_dir = Path(jobs_root) / "records"
        self.records_dir.mkdir(parents=True, exist_ok=True)

    def new(self, tool_id: str, params: dict) -> JobRecord:
        rec = JobRecord(
            uuid.uuid4().hex[:12], tool_id, params
        )  # job_id fixed at submit
        self._write(rec)
        return rec

    def _write(self, rec: JobRecord) -> None:
        (self.records_dir / f"{rec.job_id}.json").write_text(json.dumps(asdict(rec)))

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
        return [
            JobRecord(**json.loads(p.read_text()))
            for p in self.records_dir.glob("*.json")
        ]

    def reconcile_on_boot(self) -> None:
        # ADR0005: local - re-queue
        for rec in self.load_all():
            if rec.status in ACTIVE:
                self.set(rec, Status.QUEUED, handle=None)
