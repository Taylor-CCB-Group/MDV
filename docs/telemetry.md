# Telemetry from MDV

MDV can send traces, metrics and browser timings to an OpenTelemetry collector, and it sends nothing unless a deployment sets the variables below. [ADR 0025](adr/0025-self-hosted-grafana-observability.md) explains the design. The Grafana stack that receives the data runs beside MDV on each server and is kept in a separate operations repository.

## Settings

| Variable | What it does |
|---|---|
| `OTEL_EXPORTER_OTLP_ENDPOINT` | Switches server telemetry on and says where to send it, such as `http://alloy:4318` |
| `OTEL_SERVICE_NAME` | Names the deployment in Grafana. Give each deployment on a server its own name |
| `OTEL_TRACES_SAMPLER`, `OTEL_TRACES_SAMPLER_ARG` | Keep a fraction of traces, such as `parentbased_traceidratio` and `0.2` for 20%. Request metrics still count every request |
| `FARO_COLLECTOR_URL` | Switches browser telemetry on and says where MDV forwards it, such as `http://alloy:12347/collect` |
| `FARO_API_KEY` | Sent with each forwarded browser post. It must match the collector's key |
| `FARO_SAMPLE_RATE` | Records a fraction of browser sessions, such as `0.2` for 20% |
| `GEVENT_MAX_BLOCKING_TIME` | How long the event loop may stay blocked before MDV reports it, in seconds. gevent's default is 0.1 |

## Connect a deployment

Start the stack first, because it creates the `mdv-telemetry` network. Then add the following to the MDV service in the deployment's Compose file and recreate that service:

```yaml
services:
  mdv_app:
    environment:
      OTEL_EXPORTER_OTLP_ENDPOINT: http://alloy:4318
      OTEL_SERVICE_NAME: <deployment name>
      FARO_COLLECTOR_URL: http://alloy:12347/collect
      FARO_API_KEY: <the stack's FARO_API_KEY>
    networks:
      - default
      - mdv-telemetry

networks:
  mdv-telemetry:
    external: true
```

With `FARO_COLLECTOR_URL` set, MDV's pages load Grafana Faro, which records page loads, Web Vitals, request timings and uncaught errors, but not console messages. Faro posts to `telemetry/faro` on the same deployment, and MDV forwards each post to the collector. MDV refuses posts from other sites, posts over 256 KiB and more than 10 posts a second from one address, with bursts of up to 100. With auth on, sign-in is needed too.

## Names the dashboards read

The stack's dashboards and Grafana settings query these names. Renaming or removing one, or changing its labels, needs a matching change to the stack.

Server metrics, with the name Prometheus stores in brackets:

- `mdv.operation.duration` (`mdv_operation_duration_seconds`), labelled `mdv.operation` and `outcome`
- `mdv.operation.memory_growth` (`mdv_operation_memory_growth_bytes`), labelled `mdv.operation` and `outcome`
- `mdv.event_loop.blocked_time` (`mdv_event_loop_blocked_time_seconds_total`)
- `http.server.request.duration` (`http_server_request_duration_seconds`), labelled `http.route` and `http.response.status_code`, from the Flask instrumentation, which MDV switches to the current HTTP names by setting `OTEL_SEMCONV_STABILITY_OPT_IN=http`
- `process.cpu.time`, `process.memory.usage` and `process.open_file_descriptor.count`, from the system metrics instrumentation

Browser measurements, sent through Faro:

- `mdv_spatial_first_load` and `mdv_spatial_settle_load`, each with a `duration_ms` value

Logs:

- While a request is traced, every line of a log message ends with ` trace_id=<32 hex digits>`. Grafana uses it to link a log line to its trace and to find a span's logs.
