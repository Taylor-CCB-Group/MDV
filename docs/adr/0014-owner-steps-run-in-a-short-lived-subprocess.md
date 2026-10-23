# Owner steps run in a short-lived subprocess the driver waits on

**Status:** accepted (jobs framework, server integration).

## Context and decision

The driver (ADR-0012) does two pieces of heavy owner-side work for each job: it materializes the
tray from `datafile.h5`, and it ingests the worker's output columns back into the project. Both are
blocking HDF5 and numpy work in C extensions with no gevent yield point.

The production catalog runs under gunicorn's gevent worker (`-k gevent`), which calls
`monkey.patch_all()` before it loads the app. There, `threading.Thread` is a greenlet, so the driver
is a greenlet and any heavy step it runs stops every other request for as long as the step takes.
A probe that read and wrote 200 MB of gzip-compressed HDF5 stalled the hub for 4.6 s when run on a
greenlet, the full duration of the work. The same work on a native thread left the longest hub gap
at 0.02 s.

A native thread is not safe either. The driver would then open `datafile.h5` truly concurrently with
request handlers in the same process, and HDF5 refuses to open a file for writing in a process that
already has it open read-only (`OSError: ... file is already open for read-only`). `_get_h5_handle`
retries only lock errors, so the request fails. The project lock cannot prevent this, because it is
an fcntl file lock and does not exclude threads within one process. A native lock would not help,
because a greenlet blocked on it blocks the whole hub and brings the stall back.

We decided that **materialize and ingest run as owner steps, each in its own short-lived process
that the driver starts and waits on.** The driver stays where it is. Waiting on a subprocess is
cooperative under gevent, so the hub stays free, and any clash with a request becomes a
cross-process HDF5 lock error, which `_get_h5_handle` already retries. The same change fixes the
single-project and desktop servers, whose driver is a native thread and so already had the
same-process conflict.

## How an owner step is called

```
python -m mdvtools.jobs.owner_step <materialize|ingest> <project_dir> <workspace_dir>
```

- The job record is passed as JSON on stdin, so the step needs no `records_root` and the driver
  stays the only writer of job records.
- `materialize` writes the tray under the project read lock.
- `ingest` writes the output columns and their provenance pointers under the project write lock,
  then prints the provenance as JSON on stdout. The driver saves it on the DONE record.
- A non-zero exit fails the job. The last line of stderr, which is the exception's type and
  message, goes into the record's `error`; the full stderr goes to the server log.
- The driver waits for the step within its tick, so ticks stay serial. While one job stages, other
  jobs are not polled. If staging grows slow enough for that to matter, the step can become its own
  job state, launched in one tick and checked in later ones, without changing this interface.

The owner is now the web server together with these processes. The single-writer rule of ADR-0004
holds: only the owner writes the project store, and the worker never does.

## Considered options

- **Run the whole driver on a native thread** (`gevent.monkey.get_original("threading", "Thread")`).
  Frees the hub, but creates the same-process HDF5 conflict above in the catalog.
- **Keep the driver a greenlet and send the heavy steps to gevent's threadpool.** Same conflict, and
  the jobs code becomes gevent-aware.
- **Run the driver as its own long-lived process or container.** Removes both problems and is where
  similar systems end up, but submit can no longer wake the driver with an in-process event, the
  driver needs the catalog's project list, and every deployment gains a process. Owner steps are the
  step towards it: if the driver moves out, the owner steps become calls it makes directly.

## Consequences

- Each owner step pays a Python process start. That is small next to the HDF5 work it moves.
- Staging memory counts against the server's container, as before, but an out-of-memory kill now
  usually takes the owner-step process (the largest) rather than the server, so only that job fails.
- Correctness against requests relies on cross-process locks: HDF5 file locking and the fcntl
  project lock. If the project store sits on a network filesystem (NFS, Azure Files over SMB),
  check that both work there before relying on them.
- The new columns are written by another process. The server sees them because `MDVProject` reads
  `datasources.json` from disk on every access.

Prior art: Galaxy runs job handlers in processes separate from its web workers and sets dataset
metadata outside the web server; Airflow runs its scheduler apart from the webserver.
