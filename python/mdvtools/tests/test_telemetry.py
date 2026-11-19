import logging
from types import SimpleNamespace

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


SAME_SITE = {"Sec-Fetch-Site": "same-origin"}


@pytest.fixture
def browser_app(monkeypatch):
    """An app with browser telemetry on, recording what it forwards to Alloy."""
    requests = pytest.importorskip("requests")
    monkeypatch.setenv("FARO_COLLECTOR_URL", "http://alloy:12347/collect")
    monkeypatch.setenv("FARO_API_KEY", "secret")
    monkeypatch.setenv("OTEL_SERVICE_NAME", "mdv-test")
    forwarded = []

    class Upstream:
        status_code = 202

    def post(url, data, headers, timeout):
        forwarded.append((url, data, headers))
        return Upstream()

    monkeypatch.setattr(requests, "post", post)
    app = Flask(__name__)
    telemetry.register_browser_telemetry(app)
    return app.test_client(), forwarded


def test_browser_config_is_off_without_a_collector(monkeypatch):
    monkeypatch.delenv("FARO_COLLECTOR_URL", raising=False)
    app = Flask(__name__)
    telemetry.register_browser_telemetry(app)
    client = app.test_client()

    assert client.get("/telemetry/config").json == {"enabled": False}
    assert client.post("/telemetry/faro", data="{}", headers=SAME_SITE).status_code == 404


def test_browser_config_names_the_deployment(browser_app):
    client, _ = browser_app

    assert client.get("/telemetry/config").json == {
        "enabled": True,
        "app_name": "mdv-test-browser",
        "sample_rate": 1.0,
    }


def test_posts_from_this_site_reach_alloy_with_the_api_key(browser_app):
    client, forwarded = browser_app

    response = client.post("/telemetry/faro", data=b'{"meta":{}}', headers=SAME_SITE, content_type="application/json")

    assert response.status_code == 202
    [(url, data, headers)] = forwarded
    assert url == "http://alloy:12347/collect"
    assert data == b'{"meta":{}}'
    assert headers["X-API-Key"] == "secret"


@pytest.mark.parametrize(
    "headers",
    [
        {"Sec-Fetch-Site": "cross-site"},
        {"Origin": "https://elsewhere.example"},
        {},
    ],
)
def test_posts_from_other_sites_are_refused(browser_app, headers):
    client, forwarded = browser_app

    assert client.post("/telemetry/faro", data="{}", headers=headers).status_code == 403
    assert forwarded == []


def test_older_browsers_pass_when_origin_matches_the_host(browser_app):
    client, forwarded = browser_app

    response = client.post("/telemetry/faro", data="{}", headers={"Origin": "http://localhost"})

    assert response.status_code == 202
    assert len(forwarded) == 1


def test_large_posts_are_refused(browser_app):
    client, forwarded = browser_app

    response = client.post("/telemetry/faro", data=b"x" * (256 * 1024 + 1), headers=SAME_SITE)

    assert response.status_code == 413
    assert forwarded == []


def test_one_address_cannot_post_beyond_its_burst(browser_app, monkeypatch):
    client, _ = browser_app
    # Stop the clock, so no tokens refill while the posts run.
    monkeypatch.setattr(telemetry, "time", SimpleNamespace(monotonic=lambda: 1000.0))
    headers = {**SAME_SITE, "X-Forwarded-For": "203.0.113.7"}

    statuses = [client.post("/telemetry/faro", data="{}", headers=headers).status_code for _ in range(101)]
    other = client.post("/telemetry/faro", data="{}", headers={**SAME_SITE, "X-Forwarded-For": "203.0.113.8"})

    assert statuses[:100] == [202] * 100
    assert statuses[100] == 429
    assert other.status_code == 202


def test_operations_run_untouched_while_telemetry_is_off(monkeypatch):
    monkeypatch.delenv("OTEL_EXPORTER_OTLP_ENDPOINT", raising=False)

    @telemetry.traced("mdv.test")
    def double(x):
        """Doubles x."""
        return x * 2

    with telemetry.operation("mdv.test", **{"mdv.columns": 1}) as note:
        note(**{"mdv.bytes": 3})

    assert double(4) == 8
    assert double.__doc__ == "Doubles x."


class _Histogram:
    def __init__(self):
        self.records = []

    def record(self, value, attributes):
        self.records.append((value, attributes))


@pytest.fixture
def recorded(monkeypatch):
    """Telemetry on, with spans and histogram records kept in memory."""
    pytest.importorskip("opentelemetry.sdk")
    from opentelemetry.sdk.trace import TracerProvider
    from opentelemetry.sdk.trace.export import SimpleSpanProcessor
    from opentelemetry.sdk.trace.export.in_memory_span_exporter import InMemorySpanExporter

    monkeypatch.setenv("OTEL_EXPORTER_OTLP_ENDPOINT", "http://alloy:4318")
    spans = InMemorySpanExporter()
    provider = TracerProvider()
    provider.add_span_processor(SimpleSpanProcessor(spans))
    monkeypatch.setattr(telemetry, "otel_trace", SimpleNamespace(get_tracer=provider.get_tracer))
    duration, memory_growth = _Histogram(), _Histogram()
    monkeypatch.setattr(telemetry, "_operation_histograms", lambda: (duration, memory_growth))
    return spans, duration, memory_growth


def test_an_operation_records_a_span_its_duration_and_memory_growth(recorded):
    spans, duration, memory_growth = recorded

    with telemetry.operation("mdv.get_byte_data", **{"mdv.columns": 2}) as note:
        note(**{"mdv.bytes": 16})

    [span] = spans.get_finished_spans()
    assert span.name == "mdv.get_byte_data"
    assert span.attributes is not None
    assert span.attributes["mdv.columns"] == 2
    assert span.attributes["mdv.bytes"] == 16
    assert "process.memory.change" in span.attributes
    labels = {"mdv.operation": "mdv.get_byte_data", "outcome": "ok"}
    assert [attributes for _, attributes in duration.records] == [labels]
    assert [attributes for _, attributes in memory_growth.records] == [labels]
    assert memory_growth.records[0][0] >= 0


def test_a_failed_operation_is_recorded_as_an_error_and_still_raises(recorded):
    spans, duration, _ = recorded

    @telemetry.traced("mdv.upload.anndata")
    def fail():
        raise ValueError("bad file")

    with pytest.raises(ValueError):
        fail()

    [span] = spans.get_finished_spans()
    assert not span.status.is_ok
    assert duration.records[0][1] == {"mdv.operation": "mdv.upload.anndata", "outcome": "error"}


def test_a_blocking_report_names_the_mdv_code_innermost_first():
    report = [
        "=" * 80,
        "\n2026-11-12T14:29:13Z : Greenlet <Greenlet: get_data> appears to be blocked",
        "Blocked Stack (for thread id 0x1):",
        '  File "/app/.venv/lib/python3.12/site-packages/flask/app.py", line 880, in full_dispatch_request\n    rv = self.dispatch_request()\n',
        '  File "/app/python/mdvtools/server.py", line 362, in get_data\n    bytes_ = project.get_byte_data(...)\n',
        '  File "/app/python/mdvtools/mdvproject.py", line 2575, in get_byte_data\n    byte_list.append(...)\n',
        '  File "/app/.venv/lib/python3.12/site-packages/h5py/_hl/dataset.py", line 1100, in __array__\n    self.read_direct(arr)\n',
        "Info:",
        '  File "/app/python/mdvtools/other.py", line 1, in not_in_the_blocked_stack\n',
    ]

    assert telemetry._blocked_frames(report) == "mdvproject.py:2575 get_byte_data < server.py:362 get_data"


def test_each_greenlet_keeps_its_latest_request_span_until_the_next_one():
    pytest.importorskip("opentelemetry.sdk")
    greenlet = pytest.importorskip("greenlet")
    from opentelemetry.sdk.trace import TracerProvider

    provider = TracerProvider()
    telemetry._remember_request_spans(provider)
    tracer = provider.get_tracer(__name__)
    current = greenlet.getcurrent()

    with tracer.start_as_current_span("request") as request_span:
        with tracer.start_as_current_span("mdv.get_byte_data"):
            during = telemetry._request_spans.get(current)
    # gunicorn sends a file response after the request's span has ended.
    while_sending = telemetry._request_spans.get(current)
    with tracer.start_as_current_span("next request") as next_span:
        pass

    assert during == request_span.get_span_context()
    assert while_sending == request_span.get_span_context()
    assert telemetry._request_spans.get(current) == next_span.get_span_context()


def test_a_blocked_loop_warning_carries_the_blocked_requests_trace_id(caplog):
    pytest.importorskip("opentelemetry.sdk")
    from opentelemetry.trace import SpanContext

    span_context = SpanContext(trace_id=0xAB, span_id=0xCD, is_remote=False)
    original_factory = logging.getLogRecordFactory()
    try:
        telemetry._add_trace_ids_to_log_records()
        with caplog.at_level(logging.WARNING, logger=telemetry.logger.name):
            telemetry._log_blocked_loop(0.1, [], span_context)
            telemetry._log_blocked_loop(0.1, [], None)
    finally:
        logging.setLogRecordFactory(original_factory)

    assert [getattr(r, "trace_id") for r in caplog.records] == [format(0xAB, "032x"), ""]
