# One backend-owned, data-only tool spec (params + outputs): GuiSpecType-native, hydrated not serialized, with shape-routed, job-keyed idempotent outputs

**Status:** accepted — POC design decision (jobs framework). ChatMDV-ready, not
ChatMDV-integrated.

## Context & decision

"Registry-driven tools" is a POC success criterion, so each tool declares itself **once, on
the backend**, as a **serializable, data-only** spec — not a hardcoded form. The spec has
three consumers from one source of truth: the **frontend** hydrates it into live controls,
the **backend** validates submitted params against it, and a **future ChatMDV** reads it as
an LLM function-calling schema.

Full tool contract: `id` / `name` / `description` + **input contract** (needs `cells` + an
expression subgroup) + **params spec** + **output spec** + **worker entrypoint** (the registry
doubles as the dispatch table for both halves- `input_shape` -> materializer (in) and `tool_id → entrypoint` (out)).

## Params spec — GuiSpecType-native, hydrated not serialized

The spec reuses MDV's existing control vocabulary (`GuiSpecType`: `column`,
`single_category_selection`, `dropdown`, `slider`, …) with `columnType` constraints. We chose
this over JSON Schema because the **frontend consumer ships now** and gets MDV column pickers
+ live datasource binding for free; the LLM consumer is deferred, so its **`GuiSpecType →
JSON Schema` projection** (column→field-id, category→enum, slider→number) is a documented
future seam. (Outside bio, the same "one spec → form + validation + LLM" *is* JSON Schema —
MCP tool definitions are JSON Schema. We project to it later rather than adopt it now.)

This is a **new data-only spec the frontend *hydrates into* a live `GuiSpec`** — **not** a
serialization of `GuiSpec`, which is impossible: `GuiSpec` controls carry `func` callbacks and
`sourceColumn: () => FieldSpec` (functions, `charts.d.ts`). The spec carries no callbacks; the
frontend attaches them on arrival.

**Inter-param bindings: exactly two kinds, no general DSL.**
1. *options-from-another-param's-column* — `target` / `reference` choices are the categories of
   the column named by `groupby`.
2. *visibility-conditional-on-mode* — the comparison-mode selector (column-vs-value modes show
   `groupby`/`target`/`reference`; "use current selection" shows none).

The backend re-validates both at submit (are `target`/`reference` real categories of `groupby`
**in the filtered subset**; are both groups non-empty). The format may grow a third binding
when a second tool actually needs one — not before. (This is Galaxy tool-XML's
`<conditional>` / `data_column` territory, where all the complexity lives; we cap it
deliberately.)

## Safety

The registry is an **allow-list of vetted tools with typed, validated params**. A future LLM
selects a `tool_id` and fills typed blanks; it **never executes free-form code**. Backend
validation is the gate. (Standard function-calling safety — recorded for completeness, not the
novel part.)

## Output spec — shape-routed: outputs are placed by their declared axis

DGE is one job among many to come (annotation, imputation, pathway scoring, cell–cell
communication, trajectories). So the output spec is **not** a DGE-shaped contract. Each output
declares the **axis it is keyed by** — *what it is a function of* — and ingest routes it by that
shape. The framework never special-cases a job type; it only knows shapes.

| Declared shape | Output is a function of | Ingest target |
|---|---|---|
| `column(ds, cols)` | a row of an existing datasource | add column(s) to that datasource — `cells` (per-cell), `genes` (per-gene), or any other |
| `matrix(axisA, axisB, subgroups)` | an (A × B) pair | `rows_as_columns` subgroup(s) on A↔B; materialise axisB as a datasource if new |
| `new_entity(name, cols)` | a new kind of thing the job invents | a new datasource |

The rule in one line: **existing entity → columns on its datasource; a pair of axes → a matrix; a
newly-invented entity → a new datasource.** Most jobs are the first row — cell-type labels, QC,
doublet scores, HVG flags are all plain column writes (`column("cells", …)` / `column("genes",
…)`), and `concat_columns` is one too, on whichever datasource the user picked. The richer shapes
are for outputs that aren't single-valued per an existing entity.

### Worked example — `concat_columns` (the POC's first tool)

`concat_columns` joins two columns of a chosen datasource into a new text column on that same
datasource. Its output is one value per **row of that datasource** — a function of an existing
entity — so it declares a single output in the column-write shape:

```
column("cells", ["donor_tissue"])      # or column("genes", …), or any datasource
```

Ingest adds one `text` column to the named datasource. No new datasource, no matrix — the cheapest
shape in the table, and the one most future jobs (cell-type labels, QC flags, doublet scores) will
share.

**DGE (deferred).** DGE is the rich case — one value per *(gene, comparison)*, which needs the
`matrix` + `new_entity` shapes (a `contrasts` datasource + `genes × contrasts` stat matrices). That
worked example and the fuller `contrasts` design are **parked with DGE** (see ADR-0003) and are
out of POC scope. The taxonomy above already
covers them for when DGE returns.

### Output identity & idempotency (per shape)

**Output identity is the `job_id`, fixed at submit — never a timestamp generated at ingest.**
Ingest is **replace-or-skip keyed on `job_id`** (never append), applied to whichever shape the
output declared — replace the columns, replace the whole datasource, or replace the subgroup
matrix. This makes ingest **idempotent** — the property ADR-0005's recovery depends on: a
completion delivered twice (re-queue *or* at-least-once Socket.IO) yields **one** output, not
duplicates. A `job_id` also does collision-avoidance **better** than a timestamp (two submissions
in the same second collide; clock skew/timezones vanish). A human-readable label
(`dge_<groupby>_<target>_vs_<reference>`) is fine for *display*; any timestamp in it must be the
**submit time captured once in the manifest** and reused verbatim — "fixed at submit," not
`now()` at ingest.

For the **`column(ds, cols)`** shape (the POC's only one), "replace" is `remove_column` then
`add_column_to_group` — the column writer uses `create_dataset`, which errors if the dataset
already exists, so replace is remove-then-add (verified in `mdvproject.py`).

### Open (deferred with DGE) — the matrix write path is whole-matrix only

`add_rows_as_columns_subgroup` (`mdvproject.py`) writes a **complete** matrix and **refuses to
overwrite** (`if name in ds: raise ValueError`); there is no per-column append/replace. So
idempotency is cheap for **batch** matrix outputs (all comparisons from one `rank_genes_groups`
call → one matrix, replace whole) but **not yet supported for incremental** single-comparison-
per-job, which would need a subgroup replace/append. Until that exists the incremental case must
recompute the whole matrix, or fall back to `per_gene` columns / a per-run `new_entity`
datasource. **Capability gap to close before incremental matrix jobs.**

### Open (deferred with DGE) — the matrix's gene axis (OQ1)

`matrix("genes", …)` keys to the rows of a gene datasource. The `genes` datasource is the HVG
subset (`var`); DGE on `.raw` covers a *superset*, so genes outside `var` have no row to land in.
Either restrict DGE to HVG or key the matrix to a dedicated full-gene datasource. This is OQ1
(`.raw` preservation) resurfacing as a structural constraint — see ADR-0003.

## Input spec — shape-routed: trays are built by the input's declared shape

  The mirror of output-routing, on the way in. A job's inputs are staged into its workspace as a
  **tray** (ADR-0004) by an **owner-side materializer** that decodes at the boundary (`get_column`
  → real values, under a read lock) so the worker reads plain arrays and imports no MDV. The
  materializer is keyed by the input's **declared shape, not by tool** — the same discipline as
  ingest. Tools that share an input shape share a materializer; the per-tool halves are the
  **worker** (the computation) and the **spec**.

  | Declared shape | Input is | Materializer produces (the tray) |
  |---|---|---|
  | `columns(ds, cols)` | a set of columns of one datasource | one string dataset per `column` param (named by the param), scalar params as attrs |
  | `matrix(ds, layer, factors)` | an expression matrix + grouping factor(s) | the matrix + the factor arrays |

  The rule in one line: **a set of columns → one dataset per column; a matrix + factors → the
  matrix plus its grouping.** `concat_columns` is the first row — two `column` params on one
  datasource — and any later column-input tool (UMAP on chosen columns, a per-column annotation)
  **reuses the same materializer**; only its worker differs. DGE is the `matrix` row and is parked
  with DGE (ADR-0003).

  The registry dispatches the materializer by `input_shape` exactly as it dispatches the worker by
  `entrypoint` — the registry is the dispatch table for **both** halves of the courier (materialise
  in, run out). Adding a column-input tool adds a `ToolSpec` + a worker; it does **not** add a
  materializer.

## Consequences

- One source of truth eliminates the frontend/backend params drift a separate form + validator
  would invite, and makes the later LLM projection mechanical.
- **Shape-routing keeps the registry job-agnostic.** Adding a tool (annotation, imputation,
  pathway scoring, …) is *declaring its output shapes*, not writing a bespoke ingest handler.
  DGE's `contrasts` + matrix is the worked instance; the four shapes
  (`per_cell` / `per_gene` / `matrix` / `new_entity`) cover the foreseeable jobs, and a new shape
  is added only when a tool genuinely needs one.
- **Input shape-routing keeps materialisation job-agnostic too.** The owner builds trays per
    *input shape*, not per tool, so the materialiser set is a small fixed table
    (`columns` / `matrix`) shared across tools — the symmetric twin of output shape-routing. A new
    column-input tool reuses the `columns` materialiser; only a genuinely new input shape adds one.
- **Deferred — content-hash result cache.** A key of `hash(tool + params + cells + layer)`
  (Bazel's action-cache idea) enables *"you've already run this exact analysis — open the
  existing result, or run anyway?"* It is **self-invalidating** (editing cells/params changes
  the hash → clean miss) and — because we chose submission-identity — must be an **offer, not a
  silent swap** (preserve "I clicked Submit, I got a result"). Recorded in provenance now; not
  wired up in the POC.

  **`content_hash` semantics — analysis identity, not data freshness.** `content_hash` answers
  *which analysis*, never *which run* (`job_id` answers that). Two different jobs of the same
  analysis share a `content_hash` and differ in `job_id` — this is exactly why both are
  denormalized onto the output column (ADR-0009). Today the hash is **recorded but never read
  back**: a re-run is always a **new `job_id`, fresh compute**; the reuse-vs-new-job choice is the
  deferred cache's job and, by the offer-not-swap rule above, **the human's by design**. Honest
  scope of the POC hash: `content_hash(tool_id, params, input_filter_hash)` (`provenance.py`) — for
  `concat_columns` `input_filter_hash` is `None`, so it is **just the params**. It is **not** a hash
  of input *values*, so it will not notice an edit to underlying column data with unchanged params.
  It is an analysis/cache key, not a data-freshness guarantee; the future key tightens this by
  adding `cells + layer`.

Prior art: [MCP tools (JSON Schema)](https://modelcontextprotocol.io/specification/2025-06-18/server/tools),
Galaxy tool XML (`<conditional>` / `data_column`), [Stripe idempotency keys](https://docs.stripe.com/webhooks).
