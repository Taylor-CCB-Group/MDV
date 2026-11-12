import logging

import pytest
from flask import Flask

from mdvtools import telemetry
from mdvtools.logging_config import _Formatter


def _record(**extra) -> logging.LogRecord:
    record = logging.LogRecord("mdv", logging.INFO, __file__, 1, "hello", None, None)
    record.__dict__.update(extra)
    return record


def test_init_does_nothing_without_an_endpoint(monkeypatch):
    monkeypatch.delenv("OTEL_EXPORTER_OTLP_ENDPOINT", raising=False)
    app = Flask(__name__)

    assert telemetry.init_telemetry(app) is False
    # Instrumenting replaces wsgi_app on the instance and adds an after_request hook.
    assert "wsgi_app" not in vars(app)
    assert not app.after_request_funcs


def test_log_lines_carry_the_trace_id_only_when_there_is_one():
    formatter = _Formatter(fmt="%(message)s")

    assert formatter.format(_record(trace_id="ab" * 16)) == f"hello trace_id={'ab' * 16}"
    assert formatter.format(_record(trace_id="")) == "hello"
    assert formatter.format(_record()) == "hello"


def test_log_records_get_the_current_trace_id():
    pytest.importorskip("opentelemetry.sdk")
    from opentelemetry.sdk.trace import TracerProvider

    tracer = TracerProvider().get_tracer(__name__)
    original_factory = logging.getLogRecordFactory()
    try:
        telemetry._add_trace_ids_to_log_records()
        make_record = logging.getLogRecordFactory()

        with tracer.start_as_current_span("request") as span:
            inside = make_record("mdv", logging.INFO, __file__, 1, "x", None, None)
        outside = make_record("mdv", logging.INFO, __file__, 1, "x", None, None)
    finally:
        logging.setLogRecordFactory(original_factory)

    assert getattr(inside, "trace_id") == format(span.get_span_context().trace_id, "032x")
    assert getattr(outside, "trace_id") == ""
