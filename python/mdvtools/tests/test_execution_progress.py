from mdvtools.llm.execution_progress import (
    ProgressEvent,
    ProgressThrottler,
    build_heartbeat_event,
    friendly_subprocess_failure_message,
    infer_progress_event_from_output,
    parse_explicit_progress_line,
)


def test_parse_explicit_progress_line_valid_payload():
    event = parse_explicit_progress_line(
        'CHATMDV_PROGRESS:{"stage":"marker_ranking","pct":72,"msg":"Running marker ranking","step_index":2,"step_total":4}'
    )
    assert event is not None
    assert event.source == "explicit"
    assert event.stage == "marker_ranking"
    assert event.progress == 72
    assert event.step_index == 2
    assert event.step_total == 4


def test_parse_explicit_progress_line_malformed_payload():
    event = parse_explicit_progress_line('CHATMDV_PROGRESS:{"stage":}')
    assert event is None


def test_infer_progress_rank_genes_groups():
    event = infer_progress_event_from_output("scanpy.tl.rank_genes_groups(adata, 'leiden')")
    assert event is not None
    assert event.source == "inferred"
    assert event.stage == "marker_ranking"
    assert event.progress >= 70


def test_throttler_emits_stage_change_immediately():
    now = [100.0]

    def fake_clock() -> float:
        return now[0]

    throttler = ProgressThrottler(min_interval_seconds=5.0, clock=fake_clock)
    first = ProgressEvent(message="phase 1", progress=40, delta=0, stage="preprocess")
    second = ProgressEvent(message="phase 2", progress=45, delta=0, stage="neighbors")
    assert throttler.should_emit(first) is True
    assert throttler.should_emit(second) is True


def test_heartbeat_event_and_payload_keep_legacy_fields():
    event = build_heartbeat_event(12.7, stage="marker_ranking")
    payload = event.to_payload("abc-123")
    assert payload["id"] == "abc-123"
    assert payload["message"].startswith("Still computing")
    assert payload["progress"] == 0
    assert payload["delta"] == 0
    assert payload["stage"] == "marker_ranking"
    assert payload["source"] == "heartbeat"


def test_marker_ranking_progress_sequence_has_continuous_feedback():
    lines = [
        "scanpy.tl.rank_genes_groups(adata, 'clusters')",
        "still running computation...",
        "ranking complete",
    ]
    events = [infer_progress_event_from_output(line) for line in lines]
    assert events[0] is not None
    heartbeat = build_heartbeat_event(15.0, stage="marker_ranking")
    assert heartbeat.elapsed_seconds == 15


def test_friendly_subprocess_failure_message_for_sigkill():
    message = friendly_subprocess_failure_message(
        "Subprocess failed with returncode=-9. stdout_preview=''. stderr_preview=''"
    )
    assert message is not None
    assert "likely memory/resource limits" in message


def test_friendly_subprocess_failure_message_non_kill_is_none():
    message = friendly_subprocess_failure_message(
        "Subprocess failed with returncode=1. stdout_preview='x'. stderr_preview='y'"
    )
    assert message is None
