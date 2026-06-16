# Stage inputs into a per-job workspace (the courier model): sole-writer owner, relocatable worker, pluggable executor

**Status:** accepted — POC design decision (DGE jobs framework).

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
- **Single active writer by construction.** Only the owner touches the project store. Staging
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
