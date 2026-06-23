# Job workspace is ephemeral scratch outside the project store; provenance becomes durable at ingest

**Status:** accepted — POC design decision (jobs framework).

## Context & decision

The per-job workspace (ADR-0004) holds a staged input tray + intermediates + worker outputs.
It is **ephemeral scratch** under a **dedicated, configurable root**, kept **separate from the
durable project store** — one subdirectory per job keyed by `job_id`. (Industry standard:
Galaxy `job_working_directory`, Nextflow `workDir`.)

The **workspace scratch root** is **outside the project**, so `project export`/zip never
sweeps staged copies into the bundle, and so the worker can compute wherever it runs.
(Build-time note: verify the chosen root is provably outside the path the project exporter
walks.)

> **Revision (2026-06-23) — split the one root in two.** The POC originally placed a single
> `jobs_root` (records + workspaces together) as a **sibling** of the project. That breaks the
> catalog: `serve_projects_from_filesystem` / `mdv_desktop` do `os.listdir(base_dir)` and treat
> every top-level child as a project, so a sibling `jobs/` is mistaken for one. The two things
> bundled under `jobs_root` have opposite homes, so we separate them:
>
> - **Durable records** (ADR-0005) move **inside the project** → `<project>/jobs/records/`.
>   The catalog never scans a project's *subdir*, so it is no longer mistaken for a project;
>   provenance now travels with the project; HPC-safe because the owner always has its own
>   project. The exporter is taught to skip the `jobs/` subtree (`convert_to_static_page`
>   adds `jobs` to its ignore list), preserving the "never sweep job files into the bundle"
>   guarantee above.
> - **Ephemeral workspaces** stay in a **configurable scratch root outside the project**
>   (local temp by default; `$SCRATCH` / shared FS on HPC) — they must live wherever the
>   worker computes and must never be forced into the project dir.
>
> The **job_id** is the only link between a record (`<project>/jobs/records/<job_id>.json`)
> and its possibly-remote workspace (`<scratch_root>/<job_id>/`). Both are derived from the
> same id, so the owner can always correlate a durable record with its scratch.

## Scratch vs durable — what survives

- **Durable:** the **manifest / provenance** (tool / scanpy / MDV versions, params,
  `input_filter_hash`, duration, `job_id`), the **`input_filter_hash`** (pins *which* rows, for
  subset-taking tools),
  and the **ingested result**.
- **Scratch (GC'd):** the materialized tray + intermediates.

This **downgrades ADR-0004's wording**: the workspace is *not* "the reproducibility record."
Reproducibility = **manifest + `input_filter_hash` + result**. Exact re-run is possible only
while the durable store is **hash-unchanged**; the hash doubles as a **tamper-detector**
(re-derive → re-hash → match ⇒ reproduce identically; mismatch ⇒ inputs changed since the run,
exact reproduction impossible). This matches Galaxy (purges job dirs, keeps datasets +
provenance) and Nextflow (ephemeral `work/`, resume off cached *outputs*).

## Provenance durability is bound to ingest

The workspace and its manifest may live on storage **we do not control** — on HPC, cluster
scratch with an auto-purge policy. So the workspace manifest is **never** the provenance
archive. The worker writes its manifest into the workspace only because it can reach nothing
else (courier model, ADR-0004); the **owner reads it back at ingest and promotes it** —
provenance is **stamped onto the result datasource** (durable, in the project) **and** the
owner-side job log (ADR-0005). **Lineage is proven from the result, never from the scratch.**

Do this **uniformly** in local and HPC — local must not "cheat" by leaving the manifest in
`jobs-root` forever, because that breaks the one-design rule and HPC cannot do it anyway.

## Cleanup & consequences

- **POC cleanup:** clean the tray (`input/` + `work/`) on **successful** ingest; **keep on
  failure** for debugging; **no cron**. (TTL / retention sweeping is the multi-user future,
  deferred.)
- **HPC-ready:** the root can be a shared-FS path the cluster sees; heavy compute can use
  **node-local scratch** with selective stage-back (Nextflow `scratch` directive).
- **Edge — scratch purged before ingest** (e.g. a week of owner downtime vs a 3-day cluster
  purge): the outputs are gone → the job is `lost` → **re-queue** (safe, idempotent —
  ADR-0006). A failed job that never ingested still keeps the submit-time provenance from the
  owner-side record (ADR-0005).

Prior art: [Galaxy `job_working_directory`](https://docs.galaxyproject.org/en/master/admin/jobs.html),
[Nextflow work-dir / `scratch`](https://nf-co.re/docs/tutorials/storage_utilization/managing_work_directory_growth).
