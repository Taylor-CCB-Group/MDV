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
- Tempo: traces from MDV servers.

Dashboards in `grafana/dashboards/` load into the **MDV** folder when Grafana starts.

## Leave a container out

Alloy reads the logs of every container on the server. To skip one, add a drop stage to the `loki.process "containers"` block in `alloy/config.alloy`, then run `docker compose restart alloy`:

```alloy
stage.drop {
  source = "compose_project"
  value  = "<project name>"
}
```
