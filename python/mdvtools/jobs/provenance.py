import hashlib
import json
import time

import numpy as np


def content_hash(tool_id: str, params: dict, input_filter_hash: str | None) -> str:
    """Self-invalidating key: editing params or the input subset changes the hash, so a stale
    'you've already ran this' can never survive an edit. Order-independent."""
    payload = json.dumps(
        {"tool_id": tool_id, "params": params, "input_filter_hash": input_filter_hash},
        sort_keys=True,
        separators=(",", ":"),
    )
    return hashlib.sha256(payload.encode()).hexdigest()[:16]


def rows_hash(rows: np.ndarray) -> str:
    """The input_filter_hash of a selection (ADR-0005, ADR-0017): a hash of its sorted row indices,
    so the same rows give the same hash whatever order the client sent them in."""
    return hashlib.sha256(np.sort(rows).astype("<u4").tobytes()).hexdigest()[:16]


def build_provenance(rec, manifest: dict | None) -> dict:
    """Durable record of what produced an output. Identity + submit-time come from job record (fixed at submit);
    output stats from worker's manifest. Note: submitted_at is the captured submit time, never now() at ingest; completed_at
    is a plain fact, not a part of identity."""
    return {
        "job_id": rec.job_id,
        "tool_id": rec.tool_id,
        "params": rec.params,
        "input_filter_hash": rec.input_filter_hash,
        "content_hash": content_hash(rec.tool_id, rec.params, rec.input_filter_hash),
        "submitted_at": rec.created,  # submit-time; captured once
        "completed_at": time.time(),  # when ingest finished
        "output": manifest or {},  # the worker's manifest (rows, ...)
    }
