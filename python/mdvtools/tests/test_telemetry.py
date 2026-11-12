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
