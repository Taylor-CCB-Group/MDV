# MDV — Context Glossary

Canonical terms for the MDV / mdvtools packaging work. Glossary only — no implementation
detail, no decisions. Decisions that are hard to reverse live in `docs/adr/`.

## Packaging terms

- **mdvtools** — the single Python package published to PyPI (`pip install mdvtools`).
  There is one package and one version; "slim" and "full" are *consumption modes* of it,
  not separate products.

- **Slim core** — the **required** dependency set: what a bare `pip install mdvtools`
  gives you. Targeted at the majority of users (manipulating / viewing MDV projects).
  Includes the spatial feature (spatialdata is a required dependency, not an extra).

- **Extra** — an optional dependency group declared under
  `[project.optional-dependencies]`, installed on request via `mdvtools[<name>]`. Heavy /
  minority features only.

- **`app` extra** — the single extra. Pulls the entire optional cluster: database/server
  (sqlalchemy, psycopg2, gunicorn), chat/LLM (langchain, faiss, matplotlib), and auth
  (authlib, jose, auth0, redis). These three are a runtime import *cycle*
  (`auth → dbutils → llm → auth`), so they cannot be installed independently — they are
  one bundle. `mdvtools[app]` == the full app; there is no coherent middle ground between
  slim core and `app`.

- **Full install** — slim core **plus the `app` extra**. What the root `Dockerfile` builds
  and runs (database, view gallery, chat, auth). Never the default for a PyPI user.

- **Guarded import** — a top-level heavy import in a leaf module wrapped so that, on a slim
  install, calling the feature raises a friendly "install `mdvtools[<extra>]`" error
  instead of an `ImportError` at import time.

## Tooling terms

- **Build backend** — the tool that turns the source tree into a wheel/sdist
  (e.g. `uv_build`, `hatchling`, `poetry-core`). Determines how files are selected into
  the wheel.

- **Installer / lockfile manager** — the tool that resolves, locks, and installs
  dependencies into an environment (e.g. `uv`, `pip`, `poetry`). Produces the lockfile
  (`uv.lock`). Independent of the build backend.

## Running-jobs terms

Vocabulary for the async analysis-jobs work. The first tool is a trivial **concat-columns** job —
the framework is the deliverable; DGE is deferred (see ADR-0003/0006). Glossary only — the
decisions live in `docs/adr/` (0004–0007).

- **cells datasource** — the `obs` table; one row per cell.

- **genes datasource** — the `var` table; one row per gene.

- **expression matrix / subgroup** — the genes×cells values, stored as a rows-as-columns
  subgroup under the `cells` datasource. `gs` holds `adata.X`; each named layer holds an
  `adata.layers[...]`.

- **filter / subset** — the user's current selection over `cells` rows.

- **resolved selection** — the concrete list of `cells` rows currently selected, computed
  **on the client** and sent as a **job**'s subset — not a re-evaluable filter expression.

- **input filter hash** — a hash of the **resolved selection** that pins the exact subset a
  **job** ran on (the "subset hash" recorded in the **manifest**).

- **job** — one run of a **tool** against a datasource, a set of params, and — for tools that
  take one — an optional **subset**. (`concat_columns` takes no subset; it runs over the whole
  datasource.)

- **tool** — a registry-defined analysis with a params spec and an output spec. The first is
  `concat_columns` (joins two columns into a new text column); `dge_scanpy` is deferred.

- **executor** — the pluggable transport that runs a **job** (local subprocess now, HPC
  later).

- **worker** — the environment-agnostic compute process; reads and writes **only** its
  **workspace**.

- **owner** — the web server; the sole reader/writer of the project store.

- **workspace** — the per-job directory holding the materialized inputs, intermediates, and
  outputs for one **job**. Fixed layout: `input/`, `work/`, `output/`, terminal marker.

- **tray** — the materialized inputs the **owner** writes into the **workspace**'s `input/` for
  one **job**. HDF5, with a per-**tool** schema; decoded from the store so the **worker** reads
  plain arrays, never MDV internals. The format is an owner↔worker contract — substitutable.

- **ingest** — the **owner** step that reads **worker** outputs and writes them into the
  project.

- **job record** — the durable owner-side file for one **job**
  (`<project>/jobs/records/<job_id>.json`): status, params, and — once the job reaches DONE —
  its **provenance**. The single source of truth for a run (ADR-0005).

- **manifest** — the output stats the **worker** writes into its **workspace**
  (`output/manifest.json`): rows produced and similar facts the owner cannot know until compute
  finishes. Read back at **ingest** and embedded into **provenance**; never the durable archive
  itself (the scratch may be purged — ADR-0007).

- **provenance** — the full lineage of an output, built by the **owner** at **ingest** and stored
  as a field on the **job record**: `job_id`, `tool_id`, `params`, **input filter hash**,
  **content hash**, submit/complete times, and the embedded **manifest**. Stamped onto the
  **job record** AND, by reference, onto the output column (ADR-0007).
  _Avoid_: lineage, audit record.

- **content hash** — `hash(tool_id, params, input filter hash)`: the **analysis identity** — *which
  analysis*, not *which run* (that is the `job_id`). Two different **jobs** of the same analysis
  share a content hash. It keys the deferred result cache (ADR-0006); it is **not** a hash of input
  *values*, so it is an analysis/cache key, not a data-freshness guarantee.
  _Avoid_: job hash, cache key (when ambiguous with **input filter hash**).

- **provenance pointer** — the denormalized reference stamped on an **output column** in
  `datasources.json`: `kind` / `job_id` / `tool_id` / `content hash`. It is **not** a copy of the
  **provenance**; it dereferences to the **job record** by `job_id` (ADR-0009). Lets a UI label the
  column without opening the record.