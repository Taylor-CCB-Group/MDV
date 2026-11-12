"""OpenTelemetry for the MDV server, off unless a deployment asks for it.

``init_telemetry(app)`` does nothing unless ``OTEL_EXPORTER_OTLP_ENDPOINT`` is
set. When it is set, the server sends traces and metrics over OTLP/HTTP to that
address, normally the Alloy collector in ``observability/`` (ADR 0023):

- a trace for every Flask request, with a span for each database query
- process CPU, memory, open files, threads and garbage collection
- request duration metrics from the Flask instrumentation

Every other setting comes from the standard OpenTelemetry variables, such as
``OTEL_SERVICE_NAME`` to tell deployments apart, and ``OTEL_TRACES_SAMPLER``
with ``OTEL_TRACES_SAMPLER_ARG`` to keep a fraction of traces.

Logs are not sent over OTLP. They stay on stdout, where Alloy reads them, and
``logging_config`` adds ``trace_id=...`` to each line written during a request
so Grafana can link the line to its trace.
"""

import logging
import os

from flask import Flask, Response, request

from mdvtools.logging_config import get_logger

logger = get_logger(__name__)

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
_EXCLUDED_URLS = "/static/,/assets/,/favicon.ico,/socket.io/"

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
        _initialised = True

    FlaskInstrumentor().instrument_app(app, excluded_urls=_EXCLUDED_URLS)
    app.after_request(_name_span_after_project_route)
    logger.info(
        "Telemetry on: sending to %s as service %s",
        os.environ["OTEL_EXPORTER_OTLP_ENDPOINT"],
        os.getenv("OTEL_SERVICE_NAME", "unknown_service"),
    )
    return True


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
