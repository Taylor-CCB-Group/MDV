# Each server runs its own self-hosted Grafana stack for observability

Status: proposed

## Context and decision

MDV records almost nothing about how it runs. Python logs go to stdout as plain text, and Docker keeps them in `json-file` logs capped at 3 files of 10 MB, so a busy day pushes older lines out. The server has no metrics. The only browser timing is `src/react/spatialdata/perf.ts`, which someone has to switch on in their own browser and which writes to the console. To see why a server is slow, someone reads `docker logs` on the server.

The server's shape makes that worse. Gunicorn runs one gevent worker with no request timeout, so one slow CPU-bound call, such as reading columns from the h5 file or converting spatial data, holds up every user on that deployment, and a request that hangs is never ended.

We want the logs, CPU and memory use for each container and the server, timings for MDV's own code paths, and how long pages, data and images take to load in users' browsers. We do not want to pay for it, and we do not want to build a monitoring product.

We decided that MDV emits OpenTelemetry from the Python server and uses Grafana Faro, Grafana's browser monitoring library, in the frontend. Each server runs its own copy of Grafana Alloy, Loki, Prometheus, Tempo and Grafana. The stack is kept in a separate operations repository with the rest of the servers' setup, and `docs/telemetry.md` describes what MDV sends to it.

## How it works

### The stack on each server

The stack's Compose file starts it as the Compose project `mdv-observability`, on a Docker network named `mdv-telemetry`:

- Alloy collects everything and sends it on.
- Loki stores logs for 30 days.
- Prometheus stores metrics for 30 days.
- Tempo stores traces, the timing of each step inside a request, for its default of 14 days, and turns them into request rate, error and duration metrics in Prometheus.
- Grafana shows all three.

Alloy reads every container's logs through the Docker socket and labels each line with its Compose project and service. A server that runs several MDV deployments therefore shows each one separately in Grafana, and the deployments need no knowledge of each other or of the stack. Alloy's built-in cAdvisor and unix exporters give CPU, memory, disk and network figures for each container and for the server.

### The MDV server

The MDV image carries the OpenTelemetry SDK and sends nothing until a deployment sets `OTEL_EXPORTER_OTLP_ENDPOINT`. A deployment that wants telemetry joins the `mdv-telemetry` network and points that variable at Alloy. MDV then sends traces and metrics over OTLP. Its logs stay on stdout, where Alloy already reads them, and carry the trace ID so Grafana can link a log line to its trace. `pip install mdvtools` users and deployments without the variable see no change.

### The browser

Faro runs in the frontend and posts page load timings, Web Vitals, request timings and JavaScript errors to a route on the same MDV deployment. MDV forwards each post to Alloy's Faro receiver inside the Docker network, adding an API key that Alloy checks. Before forwarding, MDV refuses posts from another origin, posts over a size limit and posts from an IP address that has sent too many. When the deployment has auth on, MDV also requires sign-in. Alloy limits the rate again for each app, and it never downloads sourcemaps, because a client could otherwise make Alloy fetch any address the server can reach.

### Who can see it

Grafana and Alloy's own UI listen only on `127.0.0.1` on the server. People open Grafana through an SSH tunnel and sign in with the admin password from the stack's `.env` on that server. Anonymous access and sign-up are off. Nothing in the stack is reachable from outside the server, and the only public part is MDV's own forwarding route.

## Why

Every part is free to self-host. Grafana, Loki and Tempo are licensed AGPL-3.0, and Prometheus, Alloy and the Faro web SDK are licensed Apache-2.0. The AGPL only asks us to publish changes if we modify those programs and offer them over a network, and we run them unmodified.

OpenTelemetry is the vendor-neutral standard for instrumentation. If we later move to another backend that accepts OTLP, MDV's instrumentation stays as it is.

One stack per server keeps each dashboard to that server's services, keeps the data on the server it came from, and needs no network route between servers. People with access to a server already have SSH, so the tunnel needs no new accounts.

Forwarding browser data through MDV works on deployments served under a path without changing each server's reverse proxy. Alloy is never exposed, and MDV can apply the deployment's own sign-in rules.

The stack is server setup, so it lives with the servers' other setup, such as deploy and backup scripts, and each server's own settings are kept next to it. Passwords stay on each server in the stack's `.env`. MDV's repository keeps only what MDV emits. The dashboards read metric names that MDV defines, so renaming one needs a matching dashboard change in the operations repository. `docs/telemetry.md` lists those names, and `AGENTS.md` asks anyone changing them to say so.

## Considered and rejected

- SigNoz. It also runs from one Compose file and is free to self-host, but it stores everything in ClickHouse, and Faro and Alloy fit the Grafana stack directly.
- Storing metrics in MDV's database and charting them in MDV-Admin. That builds the monitoring product we want to avoid and puts high-volume writes on the app's database.
- Dozzle for logs and Beszel for container metrics. They need no code, but they cannot show MDV's own timings or anything from the browser.
- Grafana Cloud. Its free tier has limits, and the data would leave our servers.
- Sentry or GlitchTip for errors. Faro already reports browser errors, and traces show server errors.
- Mimir in place of Prometheus. It is built for clusters, and one Prometheus is enough for one server.
- Exposing Alloy's Faro receiver through each server's reverse proxy. Every server's proxy config would change, and anyone could post to Alloy directly.
- Sending Python logs over OTLP as well. Alloy already reads stdout, so every line would be stored twice.
- Keeping the stack in this repository. A metric rename and its dashboard would change in one pull request, but the servers' setup would be split across 2 repositories, and each server's own settings would have no tracked home.
- One central stack for every server. Data would leave each server, and the servers would need a network route to it.

## Consequences

- Alloy runs privileged, with the Docker socket and the host's `/var/run` mounted, because cAdvisor needs them. It can control the host, so the stack must only run where the people who administer the server want it.
- cAdvisor only works on a Linux Docker host. On macOS, container CPU and memory figures stay empty, while logs, host metrics, traces and browser data still work.
- Alloy reads the logs of every container on the server, including services that are not MDV. To leave one out, add a drop rule to the stack's Alloy config.
- When Alloy first starts, it reads each container's whole log history, and Loki refuses lines older than 7 days.
- Browser data from deployments without auth is anonymous, and anyone can send made-up figures within the rate limits. Treat browser dashboards as indicators, not records.
- The stack's volumes grow with traffic until retention removes old data, so each server needs disk to spare.
- The image versions in the stack's Compose file are pinned, and someone has to update them.
