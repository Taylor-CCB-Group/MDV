# Durable, owner-side job state: recover on restart, the submit↔record race, and the idempotency dependency

**Status:** accepted — POC design decision (DGE jobs framework).

## Context & decision

A long job must survive an owner reload/deploy, and an HPC job outlives the owner entirely.
Job state is therefore **durable on disk, owner-side** — a per-job record holding lifecycle
status, the **executor handle**, params, and the `input_filter_hash`. In-memory state is a
cache/index. On boot the owner scans records and reconciles anything in flight.

Lifecycle: `queued → staging → running → ingesting → done`; off-ramps `failed` / `cancelled`
/ `stale` / `lost`.

## Owner-side, not workspace-side

The durable record lives in the **owner's** storage, **not** the workspace. On HPC the
workspace is cluster scratch with an auto-purge policy we don't control, and the executor
handle (Slurm job-id) is received by the owner from the submit response — so the owner must
hold it itself to re-attach after a restart. The worker also writes a status/terminal marker
into its workspace (it is "across town" and can reach nothing else); the owner **reads that
back and promotes it** at ingest (provenance promotion: ADR-0007).

## Recover vs re-attach — local is not HPC

A *reliable handle* is only needed when the job **outlives the owner**:

- **Local (POC):** the worker subprocess dies with the owner. On boot there is no survivor,
  so any `running`/`staging` local record is **re-queued, not re-attached** — no OS PID is
  tracked. (A bare PID is an unsafe handle anyway: PID reuse, and a child reparented to init
  on the owner's death.)
- **HPC (future):** the job genuinely outlives the owner; the handle is the **Slurm job-id**,
  durable and safe by construction — re-attach via `GET /job/{id}`. `slurmrestd` is stateless
  (no session memory between requests), so only our durable record can reconnect it. This is
  where durable re-attach earns its keep.

## The submit↔record race

There is a window between `submit()` starting the work and the owner persisting the handle.
It **cannot** be closed atomically — "fork a process" and "fsync a record" are not one
operation (the classic **dual-write problem**). We therefore do **not** "persist atomically
with submit." Instead, **write-ahead intent + reconcile**:

1. write the record first in `staging`/`queued` with **no handle**;
2. submit;
3. record the returned handle, move to `running`;
4. on boot, a record stuck in the ambiguous window (intent written, no handle recorded) is
   **re-queued**.

This is the transactional-outbox shape, and exactly what Galaxy's handler does:
`__check_jobs_at_startup` resets and requeues a job that has a runner assigned but no external
id recorded.

## Idempotency dependency (load-bearing — do not remove silently)

⚠️ **Re-queue is safe *only because* ingest is idempotent** (ADR-0006: keyed on `job_id`). A
completion delivered twice — by re-queue *or* by at-least-once Socket.IO — must not create two
result datasources. **Do not make ingest non-idempotent without revisiting this recovery
strategy.** The worker's terminal marker (`done`/`failed`) in the workspace is the **primary**
completion signal, checked before any handle polling; handle-liveness is the secondary signal
(→ `stale`).

## Consequences

- The local "die-and-requeue" choice trades resumability for simplicity — acceptable precisely
  because idempotent ingest makes the re-run free of side effects.
- Cancellation is not built (lifecycle reserves `cancelled`): later, Slurm REST `DELETE
  /job/{id}` or signalling the local subprocess.

Prior art: [Galaxy job handler](https://docs.galaxyproject.org/en/master/_modules/galaxy/jobs/handler.html),
[slurmrestd](https://slurm.schedmd.com/rest.html),
[dual-write / transactional outbox](https://www.confluent.io/blog/dual-write-problem/),
durable execution (Temporal/Restate).
