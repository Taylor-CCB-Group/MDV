# Server-side job driver: a single in-process reconciler loop that advances jobs

**Status:** accepted — POC design decision (jobs framework, server integration).

## Context & decision

The framework exposes two owner-side entry points: `JobManager.submit(tool_id, params)` and
`JobManager.tick()`. In the POC a caller loops `tick()` to walk jobs to completion. Wiring the
framework into the Flask server (for the frontend jobs selector on the Vite dev server) raises the
question ADR-0008 left open: **who drives `tick()` in the server?**

The obvious cheap answer is *request-triggered* — run `tick()` inside the status-poll endpoint. We
deliberately did **not** take it. A **job** on Slurm runs for minutes or hours and no browser polls
for that long, so a request-triggered tick would only ingest a finished job while a client happens
to be watching. Completion must happen server-side regardless of any client.

This ADR records the decision: **a single background driver advances every project's jobs,
independent of any client.** It runs as an in-process daemon thread for local dev and can move to
its own process on HPC. This is the same reconciliation loop ADR-0005 already describes (`tick()` +
boot reconcile are poll-based and idempotent); this ADR says *where its loop runs* and *how submit,
recovery, and failures fit around it*.

## Why a background driver, not request-triggered

The driver gives **liveness, not latency**: it promises each job eventually advances, on no fixed
deadline. That is the property Slurm needs, where the owner must ingest a job when `sacct` reports
it finished, hours after the last client poll. Request-triggered dispatch couples "a job progresses"
to "someone is watching", which is exactly wrong for long remote work. Notifying a watching client
promptly is a *separate* concern (ADR-0008's notification seam) — polling now, Socket.IO push later —
and is not a reason to fold job progress into the request path.

## One driver, many managers

The **driver** is a single loop per server process that iterates every project's `JobManager`
through a process-wide registry (a `JobService` singleton owning the registry, the wake **event**,
and the thread). One loop makes the "exactly one driver" invariant trivially true and visible in one
place, and bounds threads at one regardless of catalog size. Ticks are serialised, so a slow ingest
in one project delays another project's tick by a cycle; for long jobs that cross-project delay is
irrelevant, and the driver thread never blocks the web server.

The invariant "exactly one driver" has three deployment expressions that leave the framework
untouched: the single-worker dev server hosts the one thread; an HPC deployment runs one driver
process on a login/service node; a Kubernetes deployment uses one replica or a leader lease. Local
is the special case (its executor holds in-memory process handles, so the driver co-hosts the work);
Slurm and K8s are the general case (durable handles, a stateless reconciler that reattaches). We
design to the general case and let Local be the degenerate one, as ADR-0005's reconcile already
does.

### Where the JobService singleton lives (implementation note)

The `JobService` singleton is a module-level global in `server.py`, reached by a plain import from
both the `POST /jobs` route and the server-boot code that runs the recovery scan and starts the
thread. This was grilled against the alternative of attaching it to the Flask app object.

In production the two are equivalent, because there is exactly one Flask app per process:
single-project builds its own app, multi-project shares one app built in `mdv_server_app.py`. The two
designs diverge only when one process builds several apps, which happens only in the test suite (each
test calls `build_app` again).

Three points settled it for the module global:

- The driver thread runs outside any Flask application or request context, so `current_app` is not
  available there. App-attachment would force the boot code and the thread to capture and pass the
  concrete app object anyway, which is the same explicit global reference the module singleton
  already is.
- Multi-project mode calls `build_app` once per project on the shared app, so app-attachment needs an
  attach-once guard (create the service only if the app does not already hold one). That guard is
  itself a singleton keyed on the app, so app-attachment does not remove global state, it relocates
  it.
- Worker placement is orthogonal. The `JobService`, the driver, and the managers stay in the server
  process regardless; only the worker (`run_worker`) runs on a Slurm compute node via the executor,
  over the shared filesystem (ADR-0010). Where the singleton lives has no bearing on remote execution.

The one cost of the module global is that tests in a single process share it. A pytest fixture that
resets the global between tests restores isolation, which is cheaper than the app-attachment plumbing
it would replace.

## Submit is write-ahead; the driver owns dispatch

`submit()` writes the `QUEUED` **job record** (write-ahead intent, ADR-0005) and sets the wake
**event**. It does **no** dispatch. The driver owns all dispatch (materialise the **tray**, call
`executor.submit`) and all **ingest**. Three reasons:

- **Keep heavy work off the hub.** The server runs on gevent (`gevent.pywsgi.WSGIServer`), a
  cooperative green-thread model. Materialising a tray is blocking HDF5 work in a C extension with no
  yield point, so doing it in the request greenlet would freeze every other client for its duration.
  The driver is a real OS thread — the codebase does not monkeypatch gevent, so `threading.Thread`
  is genuine — and h5py/numpy release the GIL during native work, so the heavy step there does not
  block the hub.
- **Single writer.** Because the driver is a real OS thread and the request is a greenlet, they run
  concurrently and the OS can interleave them. If both dispatched, they would race the JobStore and
  the executor's process handles. One dispatcher removes that race by construction.
- **Forward-compatible.** When the driver moves to its own process, `submit()` *cannot* dispatch —
  the executor lives in the driver. Deciding it this way now makes that a deployment change, not a
  rewrite.

The **event** doing double duty (a submit nudge plus the loop's idle bound) gives near-inline start
latency without putting the materialise on the request path: `submit()` nudges, the driver wakes at
once and dispatches on its own thread. The only thing given up is a submit response that literally
says `RUNNING`; the transition happens a beat later on the driver, and submit returns `202` with the
`job_id`.

## The loop

A plain daemon `threading.Thread` (dies with the process) runs one loop with two modes:

- while any manager has an active job (`QUEUED`/`STAGING`/`RUNNING`/`INGESTING`), wait on the event
  with a small configurable interval (~1s dev; widen for HPC, where `squeue`/`sacct` calls cost),
  then advance every manager. The interval governs only how fast a finished job is noticed;
  completion arrives with no event to nudge us (a subprocess exits, or the worker writes its marker).
- while nothing is active anywhere, wait on the event with no timeout, parking until a submit nudges.

This makes idle projects cost nothing and needs no per-manager "active set" bookkeeping until the
catalog is large enough for the per-tick directory listings to matter.

## Recovery scan at startup

The driver's first action is a **recovery scan**: sweep the catalog, cheaply filter each project for
active job records (one `stat` per project, skipping instantly when a project has no `jobs/`), and
build plus reconcile a manager only for the projects with in-flight work. `get_or_create` lazily
builds a manager when a live request needs one for a project the scan did not cover. This gives
boot-time reconcile — jobs reattach or re-queue after a restart with no client having to open the
project — at genuinely O(N) cheap cost, because the filter runs before the expensive manager build.

This is the standard startup-recovery pattern: `slurmctld` recovers its queue from `StateSaveLocation`
and reattaches to running jobs; a database replays its write-ahead log to finish in-flight
transactions; a Kubernetes controller LISTs existing objects and reconciles before it WATCHes. The
filter is what keeps it cheap. If N ever grows enough that even the stat sweep matters, the next step
is a single top-level index of projects-with-active-jobs, at the cost of keeping that index
consistent — deferred until profiling asks for it.

## Concurrency and locking (minimal)

`MDVProject` already carries a `fasteners.InterProcessReaderWriterLock` (`project.lock("read"|"write")`)
that works across threads and processes — the same lock covers the dev daemon thread now and a
separate driver process later. Its discipline is applied unevenly today (only `rename_view` and the
tray materialisers take it). For the driver we:

- wrap ingest's project mutation in `project.lock("write")` — the one place the driver writes the
  project.
- rely on HDF5's own file locking plus the `_get_h5_handle` retry loop (built for multi-accessor
  access) to mediate an ingest write against a concurrent `/get_data` read of the same `datafile.h5`.
- make `JobStore` record writes atomic (temp file plus `os.replace`), so the driver's `load_all` on
  its thread can never read a half-written record a submit greenlet is mid-write.

**Known gap:** `set_column` also does a read-modify-write of the datasource metadata JSON, which
neither HDF5 locking nor the write lock covers unless every other metadata writer also takes the
lock. In single-user dev the lost-update window is small. Completing the read-lock discipline across
the server's data routes is separate hardening, tracked as its own issue.

## Failure handling

The status marker is the **worker's** signal about the worker's result; `tick()` reads it and
`poll` backs it up for a worker that died without writing one (ADR-0005). An **owner-side** exception
is different: ingest's `set_column` throws, or materialise throws, though the worker may have written
`done`. There is no marker for it — but the response is the *same* JobStore transition to `FAILED`
the marker path already makes. So:

- catch owner-side exceptions **per record** inside `tick`, set that record `FAILED` with a new
  `error` field on `JobRecord`, log the traceback, and continue. The record is the source of truth
  the status endpoint reads; the log is the detail. No log scraping to learn which job failed,
  because the record is in hand at the catch point.
- errors that cannot be pinned to a single record (a whole `load_all` failing) log-and-continue at
  the manager level.

Failing on the first owner-side exception gives up a free idempotent retry (a failed ingest leaves
the workspace intact and the marker still `done`). We accept that for the POC because exceptions that
escape ingest are usually terminal, and the common transient one — h5 lock contention — is already
retried inside `_get_h5_handle`. Bounded-retry-then-`FAILED` (allow K owner-side errors before giving
up) is the hardening that reclaims the rare transient without changing the "record is the source of
truth" contract.

## Recovery-scan corruption: quarantine, don't halt and don't bury

A corrupt project found during the recovery scan is a distinct risk class, and the two easy answers
are both wrong. Halting the whole scan (a database's fail-fast on corruption) would deny recovery to
every *other* project. Logging and moving on buries a corruption event that needs a human. The
patterns that handle this — Kubernetes workqueues (per-item backoff, `Forget` after max retries,
record the failure on the object), dead-letter queues (move the poison message aside and alert), and
crash-only/microreboot (confine the failure to the smallest component) — converge on isolate-plus-
surface, split by whether the error is transient or terminal.

So the policy:

- **Parse per file, not per project.** `load_all` wraps each record's parse in its own guard, so one
  malformed file loses only itself; the project's other records still recover.
- **Quarantine the bad file**, do not delete or silently skip it — move it aside (for example a
  `records/quarantine/` folder) so it leaves the active set but is preserved for inspection.
- **Surface it durably** on a health/status signal (a count of quarantined records), not only in the
  log.
- **Split transient from terminal by exception type.** A structural/parse error is terminal:
  quarantine, stop retrying (the file will not heal itself). An IO/lock error is transient: leave it
  for the next scan.
- **Never abort the scan.** A corrupt project sits quarantined while every other project recovers.

## HTTP surface (context)

The selector talks to per-project routes registered in `add_project` (the existing closure-over-
`project` convention): `GET /project/<id>/jobs/tools` serialises the registry for the form,
`POST /project/<id>/jobs` returns `202` with `{job_id}` after re-validating (ADR-0006) and nudging,
and `GET /project/<id>/jobs` / `GET /project/<id>/jobs/<job_id>` are **pure reads** of a client-safe
record view (internal `handle` and workspace paths dropped). Status is a pure read precisely because
the driver, not the request, advances jobs.

## Consequences

- The server gains one daemon thread and a `JobService` singleton; `submit()` loses its inline
  `_dispatch()` and the driver owns dispatch, so heavy work stays off the gevent hub and there is a
  single writer of executor and project state.
- `JobRecord` gains an `error` field; `JobStore` writes become atomic. Both are backward-compatible
  (old records default `error=None`).
- Boot-time reconcile is free for projects with no jobs and cheap for the rest; a restart no longer
  needs a client to reopen a project for its jobs to advance.
- The driver moving to its own process (HPC/K8s) is a deployment change, because submit already only
  enqueues and the reconciler is already stateless for durable-handle backends.
- Two follow-ups are recorded, not built: completing the read/write lock discipline across the data
  routes, and bounded-retry-then-`FAILED` for owner-side exceptions.

Prior art:
[Kubernetes client-go workqueue](https://pkg.go.dev/k8s.io/client-go/util/workqueue),
[controller-runtime rate limiting](https://danielmangum.com/posts/controller-runtime-client-go-rate-limiting/),
[dead-letter queues / poison messages](https://www.glukhov.org/app-architecture/integration-patterns/dead-letter-queues/),
[PostgreSQL zero_damaged_pages](https://postgresqlco.nf/doc/en/param/zero_damaged_pages/),
[crash-only software](https://en.wikipedia.org/wiki/Crash-only_software),
[microreboot](https://en.wikipedia.org/wiki/Microreboot),
[Airflow executors](https://www.astronomer.io/docs/learn/airflow-executors-explained/).
