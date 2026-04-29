from __future__ import annotations

import json
import re
import time
from dataclasses import dataclass
from typing import Any


PROGRESS_PREFIX = "CHATMDV_PROGRESS:"


@dataclass(frozen=True)
class ProgressEvent:
    message: str
    progress: int
    delta: int
    stage: str | None = None
    step_index: int | None = None
    step_total: int | None = None
    eta_seconds: int | None = None
    elapsed_seconds: int | None = None
    source: str = "inferred"

    def to_payload(self, request_id: str) -> dict[str, Any]:
        payload: dict[str, Any] = {
            "message": self.message,
            "id": request_id,
            "progress": self.progress,
            "delta": self.delta,
            "source": self.source,
        }
        if self.stage is not None:
            payload["stage"] = self.stage
        if self.step_index is not None:
            payload["step_index"] = self.step_index
        if self.step_total is not None:
            payload["step_total"] = self.step_total
        if self.eta_seconds is not None:
            payload["eta_seconds"] = self.eta_seconds
        if self.elapsed_seconds is not None:
            payload["elapsed_seconds"] = self.elapsed_seconds
        return payload


_SCANPY_PATTERNS: list[tuple[re.Pattern[str], tuple[str, int, str]]] = [
    (
        re.compile(r"rank_genes_groups", re.IGNORECASE),
        (
            "marker_ranking",
            72,
            "Running marker ranking across clusters (this can take a few minutes on large datasets)",
        ),
    ),
    (re.compile(r"\bneighbors?\b", re.IGNORECASE), ("neighbors", 60, "Computing neighborhood graph")),
    (re.compile(r"\bumap\b", re.IGNORECASE), ("umap", 70, "Computing UMAP embedding")),
    (re.compile(r"highly_variable_genes", re.IGNORECASE), ("highly_variable_genes", 50, "Selecting highly variable genes")),
    (re.compile(r"normalize_total|log1p", re.IGNORECASE), ("preprocess", 45, "Preprocessing expression values")),
]


def parse_explicit_progress_line(line: str) -> ProgressEvent | None:
    if not line.startswith(PROGRESS_PREFIX):
        return None
    raw = line[len(PROGRESS_PREFIX):].strip()
    try:
        obj = json.loads(raw)
    except json.JSONDecodeError:
        return None
    msg = str(obj.get("msg") or obj.get("message") or "Progress update")
    progress = int(obj.get("pct", obj.get("progress", 0)))
    progress = max(0, min(100, progress))
    stage = obj.get("stage")
    delta = int(obj.get("delta", 0))
    event = ProgressEvent(
        message=msg,
        progress=progress,
        delta=max(0, delta),
        stage=str(stage) if stage is not None else None,
        step_index=_as_int_or_none(obj.get("step_index")),
        step_total=_as_int_or_none(obj.get("step_total")),
        eta_seconds=_as_int_or_none(obj.get("eta_seconds")),
        elapsed_seconds=_as_int_or_none(obj.get("elapsed_seconds")),
        source="explicit",
    )
    return event


def infer_progress_event_from_output(line: str) -> ProgressEvent | None:
    text = line.strip()
    if not text:
        return None
    for pattern, (stage, pct, message) in _SCANPY_PATTERNS:
        if pattern.search(text):
            return ProgressEvent(
                message=message,
                progress=pct,
                delta=0,
                stage=stage,
                source="inferred",
            )
    return None


class ProgressThrottler:
    def __init__(self, *, min_interval_seconds: float = 3.0, clock: callable = time.monotonic):
        self.min_interval_seconds = min_interval_seconds
        self.clock = clock
        self._last_emit_ts = 0.0
        self._last_stage: str | None = None

    def should_emit(self, event: ProgressEvent) -> bool:
        now = self.clock()
        stage_changed = event.stage is not None and event.stage != self._last_stage
        if stage_changed or (now - self._last_emit_ts) >= self.min_interval_seconds:
            self._last_emit_ts = now
            self._last_stage = event.stage
            return True
        return False


def build_heartbeat_event(
    elapsed_seconds: float,
    *,
    stage: str | None,
    step_index: int | None = None,
    step_total: int | None = None,
) -> ProgressEvent:
    elapsed_int = max(0, int(elapsed_seconds))
    if step_index is not None and step_total is not None:
        prefix = f"Step {step_index}/{step_total}"
    else:
        prefix = "Analysis"
    stage_text = stage or "analysis"
    return ProgressEvent(
        message=f"{prefix}: still running {stage_text} ({elapsed_int}s elapsed)",
        progress=0,
        delta=0,
        stage=stage,
        step_index=step_index,
        step_total=step_total,
        elapsed_seconds=elapsed_int,
        source="heartbeat",
    )


def _as_int_or_none(value: Any) -> int | None:
    try:
        if value is None:
            return None
        return int(value)
    except (ValueError, TypeError):
        return None


def friendly_subprocess_failure_message(stderr_diagnostic: str) -> str | None:
    """
    Convert kill-style subprocess diagnostics into a user-friendly explanation.
    Returns None for generic failures that should keep default error text.
    """
    match = re.search(r"returncode=(-?\d+)", stderr_diagnostic)
    if not match:
        return None
    code = int(match.group(1))
    if code not in (-9, 137):
        return None
    return (
        "This analysis was stopped by the runtime (likely memory/resource limits), "
        "not by a Python code error. This is common with very large datasets "
        "(for example marker-gene ranking across ~1M cells). "
        "Try a lighter run (subset/downsample/per-group) or rerun in a higher-memory environment."
    )
