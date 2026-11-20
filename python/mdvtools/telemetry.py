"""OpenTelemetry for the MDV server, off unless a deployment asks for it.

``init_telemetry(app)`` does nothing unless ``OTEL_EXPORTER_OTLP_ENDPOINT`` is
set. When it is set, the server sends traces and metrics over OTLP/HTTP to that
address, normally the Alloy collector in ``observability/`` (ADR 0025):

- a trace for every Flask request, with a span for each database query
- process CPU, memory, open files, threads and garbage collection
- request duration metrics from the Flask instrumentation

Every other setting comes from the standard OpenTelemetry variables, such as
``OTEL_SERVICE_NAME`` to tell deployments apart, and ``OTEL_TRACES_SAMPLER``
with ``OTEL_TRACES_SAMPLER_ARG`` to keep a fraction of traces.

Logs are not sent over OTLP. They stay on stdout, where Alloy reads them, and
``logging_config`` adds ``trace_id=...`` to each line written during a request
so Grafana can link the line to its trace.

``register_browser_telemetry(app)`` adds the routes that pass browser telemetry
from Grafana Faro on to Alloy. It has its own switch, ``FARO_COLLECTOR_URL``.
"""

import functools
import logging
import os
import re
import time
import weakref
from collections.abc import Callable, Iterator, Sequence
from contextlib import contextmanager
from typing import TYPE_CHECKING, ParamSpec, TypeVar

from flask import Flask, Response, request

from mdvtools.logging_config import get_logger

try:
    from opentelemetry import metrics as otel_metrics
    from opentelemetry import trace as otel_trace
except ImportError:  # a slim install of mdvtools, without the app extra
    otel_metrics = None
    otel_trace = None

if TYPE_CHECKING:
    from opentelemetry.sdk.trace import TracerProvider
    from opentelemetry.trace import SpanContext

logger = get_logger(__name__)

P = ParamSpec("P")
R = TypeVar("R")

# Host and container figures come from Alloy (cAdvisor and the unix exporter),
# so the server only reports on its own process.
_PROCESS_METRICS: dict[str, list[str] | None] = {
    "process.cpu.time": ["user", "system"],
    "process.cpu.utilization": ["user", "system"],
    "process.memory.usage": None,
    "process.memory.virtual": None,
    "process.open_file_descriptor.count": None,
    "process.thread.count": None,
    "cpython.gc.collections": None,
}

# Polled by browsers and health checks, so tracing them only adds noise.
_EXCLUDED_URLS = "/static/,/assets/,/favicon.ico,/socket.io/,/telemetry/"

_initialised = False


def telemetry_enabled() -> bool:
    return bool(os.getenv("OTEL_EXPORTER_OTLP_ENDPOINT"))


def init_telemetry(app: Flask) -> bool:
    """Instrument ``app`` and start exporting, if telemetry is configured.

    Call it before the database is initialised, so SQLAlchemy engines created
    afterwards are instrumented. Returns whether telemetry was switched on.
    """
    global _initialised
    if not telemetry_enabled():
        return False

    # Current HTTP attribute names. The older ones label request metrics with
    # the client's Host header, so the number of series grows with every name
    # the server is reached by.
    os.environ.setdefault("OTEL_SEMCONV_STABILITY_OPT_IN", "http")

    from opentelemetry import metrics, trace
    from opentelemetry.exporter.otlp.proto.http.metric_exporter import OTLPMetricExporter
    from opentelemetry.exporter.otlp.proto.http.trace_exporter import OTLPSpanExporter
    from opentelemetry.instrumentation.flask import FlaskInstrumentor
    from opentelemetry.instrumentation.sqlalchemy import SQLAlchemyInstrumentor
    from opentelemetry.instrumentation.system_metrics import SystemMetricsInstrumentor
    from opentelemetry.sdk.metrics import MeterProvider
    from opentelemetry.sdk.metrics.export import PeriodicExportingMetricReader
    from opentelemetry.sdk.resources import Resource
    from opentelemetry.sdk.trace import TracerProvider
    from opentelemetry.sdk.trace.export import BatchSpanProcessor

    # The providers are process-wide and can only be set once. Gunicorn's
    # --reload starts a new process, so this only guards repeated calls.
    if not _initialised:
        resource = Resource.create()
        tracer_provider = TracerProvider(resource=resource)
        tracer_provider.add_span_processor(BatchSpanProcessor(OTLPSpanExporter()))
        trace.set_tracer_provider(tracer_provider)

        reader = PeriodicExportingMetricReader(OTLPMetricExporter())
        metrics.set_meter_provider(MeterProvider(resource=resource, metric_readers=[reader]))

        SQLAlchemyInstrumentor().instrument()
        SystemMetricsInstrumentor(config=_PROCESS_METRICS).instrument()
        _add_trace_ids_to_log_records()
        _watch_event_loop(tracer_provider)
        _initialised = True

    FlaskInstrumentor().instrument_app(app, excluded_urls=_EXCLUDED_URLS)
    app.before_request(_remember_request_span)
    app.after_request(_name_span_after_project_route)
    logger.info(
        "Telemetry on: sending to %s as service %s",
        os.environ["OTEL_EXPORTER_OTLP_ENDPOINT"],
        os.getenv("OTEL_SERVICE_NAME", "unknown_service"),
    )
    return True


# ---------------------------------------------------------------------------
# MDV's own operations
# ---------------------------------------------------------------------------

_DURATION_BUCKETS = [0.005, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 1, 2.5, 5, 10, 30, 60, 120, 300, 600]
_MEMORY_BUCKETS = [mb * 1024 * 1024 for mb in (1, 10, 50, 100, 250, 500, 1024, 2048, 4096, 8192)]


@functools.cache
def _operation_histograms():
    assert otel_metrics is not None
    meter = otel_metrics.get_meter("mdvtools")
    duration = meter.create_histogram(
        "mdv.operation.duration",
        unit="s",
        description="How long each MDV operation took",
        explicit_bucket_boundaries_advisory=_DURATION_BUCKETS,
    )
    memory_growth = meter.create_histogram(
        "mdv.operation.memory_growth",
        unit="By",
        description="How much the process's resident memory grew during each MDV operation",
        explicit_bucket_boundaries_advisory=_MEMORY_BUCKETS,
    )
    return duration, memory_growth


def _resident_memory() -> int | None:
    try:
        import psutil
    except ImportError:
        return None
    return psutil.Process().memory_info().rss


def _ignore(**attributes) -> None:
    pass


@contextmanager
def operation(name: str, **attributes: str | int | float | bool) -> Iterator[Callable[..., None]]:
    """Record one MDV operation as a span, with its duration and memory growth.

    The duration goes to ``mdv.operation.duration`` and the growth in the
    process's resident memory to ``mdv.operation.memory_growth``, both labelled
    ``mdv.operation`` and ``outcome``. The span also carries the signed change
    in ``process.memory.change``.

    Yields a function that adds attributes to the span, for values only known
    at the end, such as the number of bytes returned. Does nothing while
    telemetry is off or without the app extra.
    """
    if otel_trace is None or not telemetry_enabled():
        yield _ignore
        return

    memory_before = _resident_memory()
    start = time.perf_counter()
    outcome = "error"
    with otel_trace.get_tracer("mdvtools").start_as_current_span(name, attributes=attributes) as span:
        try:
            yield lambda **more: span.set_attributes(more)
            outcome = "ok"
        finally:
            duration, memory_growth = _operation_histograms()
            labels = {"mdv.operation": name, "outcome": outcome}
            duration.record(time.perf_counter() - start, labels)
            memory_after = _resident_memory()
            if memory_before is not None and memory_after is not None:
                change = memory_after - memory_before
                span.set_attribute("process.memory.change", change)
                memory_growth.record(max(change, 0), labels)


def traced(name: str) -> Callable[[Callable[P, R]], Callable[P, R]]:
    """Decorator form of ``operation`` for a whole function."""

    def decorate(function: Callable[P, R]) -> Callable[P, R]:
        @functools.wraps(function)
        def wrapper(*args: P.args, **kwargs: P.kwargs) -> R:
            with operation(name):
                return function(*args, **kwargs)

        return wrapper

    return decorate


# ---------------------------------------------------------------------------
# The gevent event loop
# ---------------------------------------------------------------------------


def _watch_event_loop(tracer_provider: "TracerProvider") -> None:
    """Report whenever one request holds up the gevent event loop.

    MDV serves every request from one gevent worker, so CPU-bound work, such
    as reading columns from an h5 file, stalls every other user until it
    finishes. gevent's monitor thread notices when the loop has not switched
    for ``GEVENT_MAX_BLOCKING_TIME`` seconds (gevent's default is 0.1) and
    reports again for each further period it stays blocked. Each report
    adds that period to ``mdv.event_loop.blocked_time`` and logs one warning
    naming the MDV code that was running, with the trace ID of its request.
    """
    try:
        from gevent import monkey
    except ImportError:
        return
    if not monkey.is_module_patched("socket"):
        return

    import gevent
    from gevent.events import EventLoopBlocked
    from opentelemetry import metrics
    from zope.event import subscribers

    hub = gevent.get_hub()
    gevent.config.monitor_thread = True
    monitor = hub.start_periodic_monitoring_thread()
    # Only this hub is watched. Hubs that start later in other threads, such as
    # gevent's thread pool, would otherwise each start a monitor of their own.
    gevent.config.monitor_thread = False
    if monitor is None:
        return
    # gevent also prints every greenlet and thread to stderr for each report,
    # about 40 lines. The one-line warning below replaces it.
    setattr(monitor, "_show_blocking_report", lambda *args: None)

    _remember_request_spans(tracer_provider)
    blocked_time = metrics.get_meter("mdvtools").create_counter(
        "mdv.event_loop.blocked_time",
        unit="s",
        description="Time the gevent event loop spent blocked, in steps of GEVENT_MAX_BLOCKING_TIME",
    )

    def report(blocking_time: float, info: list[str], span_context: "SpanContext | None") -> None:
        blocked_time.add(blocking_time)
        _log_blocked_loop(blocking_time, info, span_context)

    def on_event(event) -> None:
        # This runs in gevent's monitor thread. Logging there would take
        # monkey-patched locks and start a second hub in that thread, so the
        # report waits for the main loop, which runs it once it is free again.
        # By then the greenlet may have moved on, so its span is looked up now.
        if isinstance(event, EventLoopBlocked) and event.hub is hub:
            span_context = _request_spans.get(event.greenlet)
            hub.loop.run_callback_threadsafe(report, event.blocking_time, list(event.info), span_context)

    subscribers.append(on_event)


# The span of the latest request each greenlet served
_request_spans: "weakref.WeakKeyDictionary[object, SpanContext]" = weakref.WeakKeyDictionary()


def _remember_request_spans(tracer_provider: "TracerProvider") -> None:
    from greenlet import getcurrent
    from opentelemetry.sdk.trace import SpanProcessor

    class RememberRequestSpans(SpanProcessor):
        def on_start(self, span, parent_context=None) -> None:
            # A span with no parent in this process starts a request, or other
            # work, such as a job. A browser's trace makes the parent remote.
            if span.parent is None or span.parent.is_remote:
                _request_spans[getcurrent()] = span.get_span_context()

    tracer_provider.add_span_processor(RememberRequestSpans())


def _remember_request_span() -> None:
    from greenlet import getcurrent
    from opentelemetry import trace

    span_context = trace.get_current_span().get_span_context()
    if span_context.is_valid:
        _request_spans[getcurrent()] = span_context
    else:
        _request_spans.pop(getcurrent(), None)


def _log_blocked_loop(blocking_time: float, info: Sequence[str], span_context: "SpanContext | None") -> None:
    from opentelemetry import trace

    span = trace.NonRecordingSpan(span_context) if span_context else trace.INVALID_SPAN
    with trace.use_span(span, end_on_exit=False):
        logger.warning("Event loop blocked for over %ss in %s", blocking_time, _blocked_frames(info))


_FRAME = re.compile(r'File "([^"]+)", line (\d+), in (\S+)')


def _blocked_frames(report: Sequence[str]) -> str:
    """The blocked greenlet's stack from a gevent report, innermost first.

    Keeps MDV's own frames when there are any, because the library frames
    under them rarely say which request was slow.
    """
    frames: list[tuple[str, str, str]] = []
    in_stack = False
    for line in report:
        if line.startswith("Blocked Stack"):
            in_stack = True
        elif in_stack and line.startswith("Info:"):
            break
        elif in_stack:
            frames.extend(_FRAME.findall(line))
    ours = [f for f in frames if "/mdvtools/" in f[0]] or frames[-3:]
    return " < ".join(f"{os.path.basename(path)}:{lineno} {func}" for path, lineno, func in reversed(ours[-6:])) or "unknown code"


def _name_span_after_project_route(response: Response) -> Response:
    """Rename the request's span after the project route that handled it.

    The Flask instrumentation names the span when the request starts, when
    Flask has only matched the catch-all project route. ``ProjectBlueprint``
    then sets ``request.url_rule`` to the inner rule, such as
    ``/project/<project_id>/get_data``. Tempo derives its request metrics from
    span names, so without this every project request would share one name.
    """
    from opentelemetry import trace
    from opentelemetry.semconv.attributes.http_attributes import HTTP_ROUTE

    rule = request.url_rule.rule if request.url_rule else None
    span = trace.get_current_span()
    if rule and span.is_recording():
        span.update_name(f"{request.method} {rule}")
        span.set_attribute(HTTP_ROUTE, rule)
    return response


def _add_trace_ids_to_log_records() -> None:
    """Give every log record the current trace ID, or an empty string."""
    from opentelemetry import trace

    make_record = logging.getLogRecordFactory()

    def make_record_with_trace_id(*args, **kwargs):
        record = make_record(*args, **kwargs)
        context = trace.get_current_span().get_span_context()
        record.trace_id = format(context.trace_id, "032x") if context.is_valid else ""
        return record

    logging.setLogRecordFactory(make_record_with_trace_id)


# ---------------------------------------------------------------------------
# Browser telemetry
# ---------------------------------------------------------------------------

# A Faro post batches up to 50 items, which stays well under this.
_MAX_BROWSER_PAYLOAD = 256 * 1024
# Each address may post 10 times a second, with bursts of up to 100. Many
# people can share one address behind a university or company network.
_BROWSER_POSTS_PER_SECOND = 10.0
_BROWSER_POST_BURST = 100.0


def register_browser_telemetry(app: Flask) -> None:
    """Add the routes that Faro, the browser telemetry library, talks to.

    ``FARO_COLLECTOR_URL`` switches browser telemetry on. It is the address of
    Alloy's Faro receiver inside the Docker network, such as
    ``http://alloy:12347/collect``. Without it, ``telemetry/config`` tells the
    page to stay off and ``telemetry/faro`` answers 404.

    The browser posts to ``telemetry/faro`` on its own deployment, and the
    server forwards each post to Alloy with ``FARO_API_KEY``. Before
    forwarding, it refuses posts from another site, posts over 256 KiB and
    posts from an address that has sent too many. When the deployment has auth
    on, the app-wide sign-in check covers these routes too.
    """
    import requests

    buckets: dict[str, tuple[float, float]] = {}

    @app.route("/telemetry/config")
    def browser_telemetry_config():
        if not os.getenv("FARO_COLLECTOR_URL"):
            return {"enabled": False}
        return {
            "enabled": True,
            "app_name": f"{os.getenv('OTEL_SERVICE_NAME', 'mdv')}-browser",
            "sample_rate": float(os.getenv("FARO_SAMPLE_RATE", "1")),
        }

    @app.route("/telemetry/faro", methods=["POST"])
    def browser_telemetry():
        url = os.getenv("FARO_COLLECTOR_URL")
        if not url:
            return "", 404
        if not _sent_by_this_site():
            return "", 403
        if request.content_length is None or request.content_length > _MAX_BROWSER_PAYLOAD:
            return "", 413
        if not _take_token(buckets, _client_address()):
            return "", 429

        headers = {"Content-Type": request.content_type or "application/json"}
        if request.content_encoding:
            headers["Content-Encoding"] = request.content_encoding
        api_key = os.getenv("FARO_API_KEY")
        if api_key:
            headers["X-API-Key"] = api_key
        try:
            upstream = requests.post(url, data=request.get_data(), headers=headers, timeout=5)
        except requests.RequestException as e:
            logger.warning(f"Could not forward browser telemetry: {e}")
            return "", 502
        return "", upstream.status_code


def _sent_by_this_site() -> bool:
    """Whether the post came from a page on this deployment.

    Browsers set ``Sec-Fetch-Site`` themselves, and pages cannot change it.
    For older browsers without it, compare ``Origin`` with the host the request
    was sent to. A client outside a browser can fake either header, which is
    why the rate limit exists too.
    """
    fetch_site = request.headers.get("Sec-Fetch-Site")
    if fetch_site is not None:
        return fetch_site == "same-origin"
    origin = request.headers.get("Origin", "")
    hosts = {request.host, request.headers.get("X-Forwarded-Host", "")}
    return bool(origin) and origin.split("://", 1)[-1] in hosts


def _client_address() -> str:
    # The reverse proxy appends the address it saw, so only the last entry is
    # trusted. Earlier entries come from the client and can be anything.
    forwarded = request.headers.get("X-Forwarded-For", "")
    return forwarded.split(",")[-1].strip() or request.remote_addr or ""


def _take_token(buckets: dict[str, tuple[float, float]], address: str) -> bool:
    """Token bucket for each address.

    It needs no lock, because MDV runs one worker process and gevent only
    switches between requests at I/O.
    """
    now = time.monotonic()
    if len(buckets) > 10_000:
        # An address idle long enough to refill completely has no state worth keeping.
        refill_time = _BROWSER_POST_BURST / _BROWSER_POSTS_PER_SECOND
        for key in [k for k, (_, last) in buckets.items() if now - last > refill_time]:
            del buckets[key]
    tokens, last = buckets.get(address, (_BROWSER_POST_BURST, now))
    tokens = min(_BROWSER_POST_BURST, tokens + (now - last) * _BROWSER_POSTS_PER_SECOND)
    allowed = tokens >= 1
    buckets[address] = (tokens - 1 if allowed else tokens, now)
    return allowed
