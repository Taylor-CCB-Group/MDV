# MDV observability stack

This folder runs Grafana with Loki, Prometheus, Tempo and Alloy beside MDV on one server. It collects every container's logs, CPU and memory figures for each container and for the server, and the traces, metrics and browser timings that MDV sends. [ADR 0023](../docs/adr/0023-self-hosted-grafana-observability.md) explains why it is built this way.

| Service | What it does | Kept for |
|---|---|---|
| Alloy | collects logs, metrics, traces and browser data, and sends them on | |
| Loki | stores logs, including browser data from Faro | 30 days |
| Prometheus | stores metrics | 30 days (`PROMETHEUS_RETENTION`) |
| Tempo | stores traces and turns them into request metrics | 14 days |
| Grafana | dashboards and search | |

## Start the stack

You need Docker with Compose. Container CPU and memory figures need a Linux Docker host. On macOS, everything else works.

1. Copy `example.env` to `.env` in this folder and set `GRAFANA_ADMIN_PASSWORD`.
2. From this folder, run `docker compose up -d`.
3. Open http://127.0.0.1:3000 and sign in as `admin` with that password.

To stop it, run `docker compose down`. Add `-v` to delete the stored data as well.

## Connect an MDV deployment

Container logs and CPU and memory figures need nothing from MDV. For MDV's traces and request metrics, start this stack first, because it creates the `mdv-telemetry` network. Then add the following to the MDV service in the deployment's Compose file and recreate that service:

```yaml
services:
  mdv_app:
    environment:
      OTEL_EXPORTER_OTLP_ENDPOINT: http://alloy:4318
      OTEL_SERVICE_NAME: <deployment name>
      FARO_COLLECTOR_URL: http://alloy:12347/collect
      FARO_API_KEY: <the FARO_API_KEY from this folder's .env>
    networks:
      - default
      - mdv-telemetry

networks:
  mdv-telemetry:
    external: true
```

Give each deployment on a server its own `OTEL_SERVICE_NAME`, because Grafana uses it to tell deployments apart. Without `OTEL_EXPORTER_OTLP_ENDPOINT`, MDV sends nothing.

MDV records every request by default. To keep a fraction of traces on a busy server, also set `OTEL_TRACES_SAMPLER=parentbased_traceidratio` and `OTEL_TRACES_SAMPLER_ARG=0.2` for 20%. Request duration metrics still count every request.

`FARO_COLLECTOR_URL` switches on browser telemetry. MDV's pages then load Grafana Faro, which records page loads, Web Vitals, how long each request takes and uncaught errors. It does not record console messages. Faro posts to `telemetry/faro` on the same deployment, and MDV forwards each post to Alloy. MDV refuses posts from other sites, posts over 256 KiB and more than 10 posts a second from one address, with bursts of up to 100. With auth on, sign-in is needed too. Faro also adds a `traceparent` header to MDV's own requests, so a request's browser span and server spans share one trace. To record a fraction of browser sessions, set `FARO_SAMPLE_RATE`, such as `0.2` for 20%.

## Open Grafana on a server

Grafana listens only on `127.0.0.1` on the server, so open an SSH tunnel from your own machine:

```sh
ssh -L 3000:127.0.0.1:3000 <server>
```

Then open http://127.0.0.1:3000 on your machine. If port 3000 is busy on either side, set `GRAFANA_PORT` in `.env`, or tunnel from another local port, such as `-L 3001:127.0.0.1:3000`.

## Find things in Grafana

Use **Explore** in Grafana's menu:

- Loki: container logs, labelled `compose_project`, `compose_service` and `container`. Browser data from Faro carries `source="faro"`, `app_name` and `kind`, such as `{source="faro", kind="exception"}`.
- Prometheus: container figures come from cAdvisor, such as `container_memory_working_set_bytes`, labelled `container_label_com_docker_compose_project` and `container_label_com_docker_compose_service`. Server figures come from the unix exporter, such as `node_cpu_seconds_total` and `node_filesystem_avail_bytes`.
- Prometheus, from MDV: request durations in `http_server_request_duration_seconds`, labelled `job` (the deployment's `OTEL_SERVICE_NAME`), `http_route`, `http_request_method` and `http_response_status_code`. Process figures include `process_memory_usage_bytes` and `process_cpu_time_seconds_total`. Tempo adds `traces_spanmetrics_calls_total` and `traces_spanmetrics_latency_bucket`, labelled `service` and `span_name`.
- Tempo: one trace for each request to MDV, with a span for each database query. Log lines written during a request end in `trace_id=...`, and Grafana links them to their trace.

Dashboards in `grafana/dashboards/` load into the **MDV** folder when Grafana starts.

## Leave a container out

Alloy reads the logs of every container on the server. To skip one, add a drop stage to the `loki.process "containers"` block in `alloy/config.alloy`, then run `docker compose restart alloy`:

```alloy
stage.drop {
  source = "compose_project"
  value  = "<project name>"
}
```
