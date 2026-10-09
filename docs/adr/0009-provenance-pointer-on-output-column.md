# Provenance on an output column is a pointer to the job record, not an inline copy

**Status:** accepted — POC design decision (jobs framework). Completes the datasource-side half
of ADR-0007.

## Context

ADR-0007 ("Provenance durability is bound to ingest") already decided that provenance is
"stamped onto the **result datasource** … AND the owner-side job log … Lineage is proven from the
result." The POC built only the **job-log half**: at ingest the owner promotes the worker's
**manifest** into full **provenance** and stores it as a field on the **job record**
(`<project>/jobs/records/<job_id>.json`). The **datasource half** — proving lineage *from the
output column itself* — was never built. This ADR records how we build it.

## Decision — a provenance *pointer*, not an inline copy

The output column carries a **provenance pointer**: one `provenance` key on the column's metadata
dict in `datasources.json`, holding a small **reference** to the producing job — never a copy of
the full provenance.

```json
"provenance": {
  "kind": "job",
  "job_id": "395b869efa5c",
  "tool_id": "concat_columns",
  "content_hash": "b4c86567cf07cf00"
}
```

We chose a pointer over an inline copy because the **job record is the single source of truth**
for a run (ADR-0005), and `datasources.json` must not be bloated with experimental, duplicated
state that can drift from the record. The pointer dereferences to the record by `job_id`.

Each field earns its place:

- **`kind: "job"`** — discriminator. Leaves room for `"manual"` / `"import"` origins later without
  a schema change.
- **`job_id`** — the dereference key → `<project>/jobs/records/<job_id>.json` (the full
  provenance). The *only* link between column and record, mirroring how `job_id` is the only link
  between a record and its scratch (ADR-0007).
- **`tool_id`** — denormalized so a UI can **label** the column ("produced by `concat_columns`")
  without opening the record, and so the label survives the record being purged (see read seam).
- **`content_hash`** — the **analysis identity** (see ADR-0006): *which analysis*, not *which run*.
  Denormalized so a future cache / staleness check can compare columns without dereferencing.

All four are computable from the `JobRecord` at ingest, before any worker manifest exists.

**One hash, not two.** The pointer stamps the *same* `content_hash` that `provenance.py`'s
`content_hash(tool_id, params, input_filter_hash)` already computes for the record — we do **not**
add a second hashing function. Column and record therefore can never disagree, and both move
together if the cache key later grows (ADR-0006's deferred `+ cells + layer`).

## Write seam — ingester stays data-only, manager stamps

The **ingester** stays a pure data-writer: it writes the worker's result column and now *reports
what it wrote*, returning `{"manifest": …, "outputs": [(datasource, column)]}`. It does **not**
build the pointer — it only ever receives `params`, not the `JobRecord`, so it lacks `job_id` and
`content_hash` by construction. The **manager** owns the job record and provenance, so the manager
stamps the pointer, once per reported output, via the existing
`project.set_column_metadata(ds, col, "provenance", pointer)` (`mdvproject.py`).

**Commit-last ordering.** The pointer is stamped **before** the record is flipped to `DONE` — the
`DONE` transition is the single commit point:

```
set INGESTING
  ingest column data        (idempotent: remove-then-add, ADR-0006)
  build_provenance
  stamp pointer on column    (idempotent: set_column_metadata updates-or-appends)
set DONE + provenance        ← single commit
rmtree workspace             (cleanup, idempotent)
```

This keeps the pointer stamp **inside** the idempotent ingest, so it inherits ADR-0006's
replace-or-skip idempotency *and* ADR-0005's recovery: any crash before `DONE` leaves the record in
an **ACTIVE** status (`INGESTING`), which recovery re-queues and re-ingests, rewriting all three
writes. The rejected order — stamp *after* `DONE` — has a permanent failure window: a crash between
"set `DONE`" and "stamp pointer" leaves a terminal (never-replayed) record whose column has no
pointer, forever. The only cost of commit-last is a sub-millisecond transient where a *concurrent*
reader can see the record present but its provenance field not yet filled; recovery makes it
self-heal.

## Read seam — `project.get_column_provenance(ds, col)`, three states

Production always keeps records at `<project>/jobs/records/` (the `records_root` override exists
only in tests), so the resolver lives on the **project** and resolves by that convention — the
project owns both the column metadata and its own `jobs/records/` dir; no `JobManager` instance is
needed to read lineage.

The contract distinguishes **three** states, because a UI should render them differently:

| Column state | Returns |
|---|---|
| no `provenance` key | `None` — column was not job-produced (imported / manual) |
| pointer + record present | the full record, marked resolved |
| pointer + record purged | the **pointer dict** alone, marked unresolved (dangling) |

The dangling case is precisely *why* `tool_id` / `content_hash` are denormalized: when the record
is gone, the column still self-describes ("produced by `concat_columns`, full details purged")
rather than collapsing to an indistinguishable `None`. (A caller that only wants a cheap label can
read `get_column_metadata(ds, col)["provenance"]` directly — the pointer is always right there.)

## Consequences

- **Re-run = last writer.** Re-ingesting the same output column replaces the column **and** its
  pointer → the column points to the newest `job_id`. The older job's record still exists in
  `records/`, so full history is preserved as separate files. The older record now claims an output
  it no longer owns — that is fine: **records are append-only history; the column pointer is the
  current truth.**
- **Purge → dangling, not broken.** A GC'd record leaves a dangling pointer; the resolver returns
  the pointer (unresolved) and the column stays fully usable.
- **Static export.** `convert_to_static_page` excludes `jobs/` (ADR-0007), so in a static bundle
  every pointer dangles. Accepted for the POC; a future option is to inline a one-line provenance
  summary at export time.
- **Two keys, never conflated** (see ADR-0006): `job_id` = *which run* (submission identity);
  `content_hash` = *which analysis* (cache identity). Both live on the column for this reason.

Prior art: [Stripe idempotency keys](https://docs.stripe.com/webhooks) (idempotent replay over
transactionality); Galaxy / Nextflow keep provenance with the dataset, not the work dir.
