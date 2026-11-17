# Stage inputs into a per-job workspace (the courier model): sole-writer owner, relocatable worker, pluggable executor

**Status:** accepted — POC design decision (DGE jobs framework).

**Amended** by ADR-0023 (proposed, Taylor-CCB-Group/MDV#665): see *Amendment: Icechunk project storage* at the end.

## Context & decision

A job's input data must be delivered to wherever the compute runs — a local subprocess
today, an HPC node tomorrow. The web server (**owner**) **materializes** a self-contained
copy of exactly the needed slice (ADR-0003) into a **per-job workspace**; the compute
**worker** reads only its workspace, runs scanpy, and writes outputs back into the workspace;
the owner then **ingests** outputs into the project. The worker behaves as if it is always
"across town" — the **courier model**.

The worker always reads a **materialized tray** in its workspace. The word is
**"materialize," not "copy":** the tray is *derived* from the store, not a duplicate of a file
that already exists.

**The framework fixes the workspace layout, not the tray's format.** Layout is uniform — `input/`
(the tray the owner writes), `work/` (worker scratch), `output/` (worker results), plus the
terminal marker (ADR-0007). *What* goes in `input/` and *how it is serialized* is the **tool's
contract**: the owner's materialize step and the worker's read step are two halves of one tool
definition, so they agree by construction. That format never leaves the owner↔worker boundary, so
it is **substitutable** without touching the framework.

**The tray is HDF5** — the format MDV's store (`datafile.h5`) and scanpy (`.h5ad`) both already
use, so no new dependency. DGE's tray is a scanpy-loadable AnnData (a derived slice: selected cells
× all genes — the store is gene-column CSC, scanpy wants cells×genes, so there is nothing to
duplicate). `concat_columns`' tray is a plain `.h5` of the two requested columns. **The owner
decodes at the boundary:** MDV stores text columns as `uint16` codes + a manifest `values` lookup,
so the owner resolves them to real values when writing the tray — the worker reads plain typed
arrays and never imports MDV internals (the rejected "worker reads the store directly" option
below).

## Why

- **One design serves local and HPC.** The worker is uniform and identical everywhere;
  "local" is just a very short courier walk — test it on a laptop, trust it on a cluster.
  Validated in production by Galaxy/Pulsar, which stage files to a remote host with **no
  shared filesystem**.
- **Single active writer by construction.** Only the owner touches the project store: the web
  server and the short-lived processes it starts to materialize and ingest (ADR-0014). Staging
  reads under the project **read-lock** (`MDVProject.lock`), so the tray is a torn-read-free
  **snapshot** frozen at submit time — a long job can never see a half-edited matrix. (This
  snapshot isolation is *why* materializing locally is not wasteful even though the data is
  "right there.")
- **Subprocess over thread.** Threads don't generalize to HPC; "a process reading a
  workspace" does. Airflow's `LocalExecutor` likewise runs subprocesses, not threads.

## Pluggable executor — and staging is a *separate* seam

The **executor is an abstraction** (`submit` / `poll` / `locate-result`):
`LocalSubprocessExecutor` now, a Slurm/HPC executor later. (Mirrors Airflow's `BaseExecutor`
→ Local/Celery/Kubernetes; CI runners' shell/docker/kubernetes executors.)

*Where compute runs* (subprocess vs Slurm) and *how inputs arrive* are **two seams**, kept
separate by the closest prior art: Pulsar's staging action is `transfer` / `copy` / **`none`**
(`none` on a shared FS); Bazel Remote Execution stages via content-addressed CAS with dedup.
The real discriminator is **shared-filesystem-or-not**, not local-vs-remote — an HPC node with
shared scratch also wants the no-copy path; only a *true* no-shared-FS remote must transfer
bytes.

## Worker entrypoint

The executor launches the worker as a plain subprocess:

```
python -m mdvtools.jobs.run_worker "module:fn" <workspace>
```

`run_worker` is a thin shim (no MDV imports): it splits `"module:fn"`, imports the module,
and calls `fn(<workspace>)` with the per-job workspace path as its only argument. That `fn`
is the tool's worker half (e.g. `mdvtools.jobs.workers.concat_worker:run`) — it reads the
materialized tray from `<workspace>/input/`, writes results + a manifest into
`<workspace>/output/`, and writes the terminal `STATUS` marker (ADR-0007) **last**. The
contract is the workspace path and nothing else, which is exactly why the same worker runs
unchanged under a Slurm/HPC executor: swap *how* the process is launched, not *what* it is.

## Considered and rejected (for the POC)

- **Per-strategy staging now (`symlink`/`none`/in-place locally).** Buys local speed but
  splits the worker into environment-aware paths, breaking the "reads only its workspace,
  identical everywhere" property that is the whole point. Rejected for the POC; recorded as a
  **future seam** — a shared-FS/HPC deployment can add a `symlink`/`none`/content-addressed
  strategy **without touching the executor contract**. (It only ever applies to a *future*
  tool that wants a whole file verbatim; for DGE the tray is computed, not duplicated.)
- **Worker reads the live project store directly (skip the tray) when local.** Forces the
  worker to understand MDV's internal CSC layout, and forfeits snapshot isolation. Rejected.

## Consequences

- The owner↔worker boundary **is** the workspace. The worker never imports MDV internals and
  never touches the project store.
- Staged inputs are **scratch**, owned by the per-job workspace; their lifecycle and the
  durable/scratch split live in ADR-0007.
- `submit` returns a **handle** the owner persists durably (ADR-0005).

Prior art: [Pulsar staging](https://pulsar.readthedocs.io/en/latest/galaxy_conf.html),
[Airflow executors](https://www.astronomer.io/docs/learn/airflow-executors-explained/),
[Bazel Remote Execution API](https://github.com/bazelbuild/remote-apis).

## Amendment: Icechunk project storage (ADR-0023, proposed)

ADR-0023 ([#665](https://github.com/Taylor-CCB-Group/MDV/pull/665)) retires `datafile.h5` for active
projects: bulk data stays in read-only source stores, and everything MDV writes goes to the
project's Icechunk repository. For this ADR:

- **The tray format does not have to change.** The tray belongs to the tool and is scratch, never
  project storage, so `.h5ad` and `.h5` trays stay valid. One of the two reasons for HDF5 goes:
  it is no longer the format of MDV's store. The other stays: scanpy reads `.h5ad`. A tool may
  choose zarr instead without touching the framework.
- **Snapshot isolation comes from a pinned snapshot, not the read lock.** `submit` records the
  repository snapshot that is current at submit time on the job record, and materialize reads that
  snapshot through a read-only session. Source stores are immutable and need no pin. The tray is
  then a snapshot of the project *as it was at submit*, which the read lock never actually gave:
  it froze the data at materialize time. The materialize description in ADR-0006 ("under a read
  lock") reads the same way.
- **Decoding at the boundary stays.** Under ADR-0024 categorical columns are integer codes plus a
  categories array in the store instead of `uint16` codes plus a `values` list in
  `datasources.json`. The owner still resolves them to real values, and the worker still imports
  nothing from MDV.
- **"Worker reads the live store directly" stays rejected,** though one of its reasons changes. A
  worker could open a pinned snapshot read-only, so snapshot isolation no longer argues against it.
  The worker would still have to understand MDV's layout, and a remote worker with no shared
  filesystem still cannot reach the store.
