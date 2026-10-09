# Dispatch without a message broker: peer-backend executors, executor-owned concurrency, notifications as a separate seam

**Status:** accepted — POC design decision (jobs framework).

## Context & decision

The framework must dispatch a submitted **job** to wherever it runs — a local subprocess now,
an HPC scheduler (Slurm) or Kubernetes later. ADR-0004 already defines the **executor** as a thin
seam (`submit` / `poll` / `locate-result`) with `LocalSubprocessExecutor` first. This ADR records
the dispatch *policy* sitting on that seam: **no message broker, concurrency owned by the
executor, and notifications kept as a separate future seam.** The obvious reach for "many users,
many concurrent jobs" is a broker-backed task queue (Celery on RabbitMQ/Redis); we deliberately
did **not** take it, and this records why so the question doesn't get re-opened by default.

## No message broker

We do not adopt Celery / RabbitMQ / Redis for *running* jobs. A broker-backed task queue is the
right tool when the web app **is** the compute cluster — it owns a pool of worker processes and
feeds them tasks. That is not this system's future:

- **The heavy compute goes to a cluster.** On HPC the scheduler is **Slurm**; on Kubernetes it is
  the K8s control plane. Those schedulers already do the hard multi-user part — fair-share across
  users, resource-aware placement (RAM/GPU), isolation, quotas — far better than a task queue. A
  broker in front of them is a queue in front of a queue; the scheduler re-queues the work anyway.
- **The worker is sealed (courier model, ADR-0004):** it reads and writes only its workspace and
  holds no outbound connections. On a real HPC node it often *cannot* reach an external broker. A
  broker design assumes the worker is a long-lived consumer wired to the broker — the opposite of
  the courier model.
- **Durable job state is already owner-side on disk (ADR-0005),** which provides the durability a
  broker's result backend would. Leaning on a broker's state would duplicate it and partly
  contradict the "owner-side record is the source of truth" rule.

The one future where a broker *is* the right default: a single big server that runs every job as a
local process for many users, with no cluster. If that becomes the deployment, the answer is a
`CeleryExecutor` (below) — a new backend, not a rewrite.

## Concurrency is the executor's job

Bare subprocesses have no concurrency limit, so `LocalSubprocessExecutor` bounds it with a pool /
semaphore (`max_concurrent_jobs`, configurable). This stops a burst of submissions from spawning
unbounded heavy processes. On a cluster backend, concurrency is the scheduler's responsibility, not
ours — that executor just submits and lets Slurm/K8s queue.

## Executors are peer backends behind one interface

Pluggability lives in the executor *interface*, not in any one technology. Future backends are
**peer classes** implementing the same `submit` / `poll` / `locate-result` contract:

- `LocalSubprocessExecutor` — subprocess + bounded pool; handle = token; poll = terminal marker (ADR-0005).
- `SlurmExecutor` — `slurmrestd` submit; handle = Slurm job-id; poll = `GET /job/{id}`.
- `K8sExecutor` — create a Job object; handle = namespace/name; poll = watch status.
- `CeleryExecutor` — `.delay()`; handle = task-id; poll = `AsyncResult` — only if a single-big-server deployment ever wants it.

Adding a backend is writing one class, not changing the framework. This mirrors Airflow (Local /
Celery / Kubernetes executors are siblings behind `BaseExecutor`) and the bioinformatics workflow
engines (Nextflow, Snakemake, Cromwell abstract the scheduler directly — local / slurm / k8s /
cloud — none built on a broker). Note that Celery and Kubernetes are *alternatives* in that model,
not a stack: you do not reach K8s *through* Celery. So "adopt Celery to stay open to K8s" does not
hold — the interface is what stays open.

## Notifications are a separate seam

"Run the job" and "tell interested parties it finished" are different concerns. Completion today is
the worker's **terminal marker** in the workspace (the durable signal, ADR-0005) plus **Socket.IO**
pushing to the user's browser — no broker needed at a single owner instance. A message bus earns a
slot only when there are **multiple owner instances** behind a load balancer and a completion event
must reach the instance holding a given user's WebSocket. That is a notification-fan-out problem,
not a job-running one; the lighter first answer is **Redis pub/sub** (Socket.IO ships an adapter for
exactly this), with RabbitMQ reserved for richer routing or external subscribers. Out of POC scope —
recorded here so the broker question lands on the *notification* seam, never on the dispatch path.

## Consequences

- The POC ships **one** executor (`LocalSubprocessExecutor`) and **no** broker — minimum
  infrastructure, every architectural seam still exercised.
- "Don't get locked into one backend" is satisfied by the *interface*, not by adopting a queue
  early; Slurm / K8s / Celery are each a later drop-in behind it.
- If multi-instance notifications arrive, add an event-bus seam (Redis first), kept separate from
  the executor.

Prior art: [Airflow executors](https://www.astronomer.io/docs/learn/airflow-executors-explained/),
[Nextflow executors](https://www.nextflow.io/docs/latest/executor.html),
[Snakemake executor plugins](https://snakemake.github.io/snakemake-plugin-catalog/),
[Celery](https://docs.celeryq.dev/),
[Socket.IO Redis adapter](https://socket.io/docs/v4/redis-adapter/).
