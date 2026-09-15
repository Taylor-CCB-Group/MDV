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
  later). It starts the work and reports completion; it moves **no data** — that is the
  **data-movement seam** (ADR-0010).

- **shared-filesystem precondition** — a remote **executor** (e.g. Slurm) assumes the **owner**
  and the **worker**'s compute node see the same POSIX filesystem at a **matching path**, so the
  **workspace** handoff is zero-copy. Satisfied by deployment (same node, NFS/Lustre mount, managed
  cloud FS), not by code (ADR-0010).

- **data-movement seam** — the (still unbuilt) layer that gets the **tray** to the **worker** and
  outputs back when there is no shared filesystem: a no-op under the **shared-filesystem
  precondition**, a copy/stage or content-addressed transfer otherwise. Kept **separate** from the
  **executor** (ADR-0004, ADR-0010). _Avoid_: conflating with **executor**.

- **worker** — the environment-agnostic compute process; reads and writes **only** its
  **workspace**.

- **owner** — the web server; the sole reader/writer of the project store.

- **driver** — the single background loop (one per server process) that repeatedly advances every
  project's **jobs**, dispatching queued ones and ingesting finished ones, so they progress whether
  or not any client is watching. Liveness, not latency: it promises each advance eventually happens,
  on no fixed deadline. Local dev runs it as a daemon thread inside the **owner**; HPC can host it in
  its own process. Exactly one runs at a time (ADR-0012).

- **recovery scan** — the **driver**'s first action at startup: sweep the catalog, cheaply filter
  each project for active **job records**, and build a manager only for those, so in-flight **jobs**
  reattach or re-queue after a restart with no client having to open the project. Filter first, so
  the cost is one stat per project, not a manager per project. The startup half of the same reconcile
  (ADR-0005) a per-project manager runs when it is first created.

- **quarantine** — where the **recovery scan** moves a **job record** it cannot parse: aside and
  preserved for inspection (not deleted, not silently skipped), surfaced on a health signal. A
  corrupt record loses only itself; the project's other records and every other project still
  recover (ADR-0012).

- **driver** — the single background loop (one per server process) that repeatedly advances every
  project's **jobs** so they progress whether or not any client is watching. It provides liveness,
  not latency: it promises each advance eventually happens, on no fixed deadline. Local dev runs it
  as a daemon thread inside the **owner**; HPC can host it in its own process. Exactly one runs at a
  time.

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
## SpatialData.js integration terms

Canonical rendering vocabulary lives in [SpatialData.js CONTEXT.md](https://github.com/Taylor-CCB-Group/SpatialData.js/blob/main/CONTEXT.md). MDV chart code uses the same terms:

- **Render Stack** — saved draw order on chart config (`renderStack.entries`), not parallel `layerOrder` arrays
- **Stack Entry** — one ordered item (`kind: spatial | host | group`)
- **Host Overlay** — MDV deck layer referenced by `hostLayerId` in the stack, resolved at runtime
- **Runtime Attachment** — `hostLayerResolver`, tooltip hooks, `deckProps`; not serialized in `entry.props`
- **Image Layer Registry** — chart-owned runtime attachment (`SpatialDataMdvReact`): callbacks that expose loaded Viv image data (`getImageLoadedDataByElementKey`, load state) from the spatial renderer. The viewer populates it; the layer dialog consumes it via `ImageLayerContextProvider`. Bridges separate React trees (chart vs dialog portal) without duplicating image load paths.
- **Spatial Image Panel Context** — React context on [`ImageLayerPanel`](src/react/components/spatialLayers/ImageLayerPanel.tsx): wraps `ImageLayerContextProvider` (upstream loaded-image defaults) and exposes `useLayerChannelState` write API (`setChannels`, `addChannel`, `removeChannel`) to channel UI. Persisted field edits go through the hook, not zustand `channelsStore`.
- **Image Layer Runtime Bridge** — pure helpers in `image_layer_runtime.ts`: `channelId`-keyed stats cache (`domains`, `raster`), viewer parallel-array sync (`channelOptions`, loading flags), and tone → `vivLayerProps` helpers. Not a second MobX persistence bridge.
- **MobX Control Island** — layer dialog UI that patches `config.renderStack` directly
- **Layer Channel Config** — serializable image channel state on `renderStack.entries[].props.channels` (`ChannelConfig` in `@spatialdata/vis`): colors, contrast limits, visibility, selections, and extension-related data fields. Tone (brightness/contrast) is **not** part of `channels`; it persists in `entry.props.vivLayerProps`. Distinct from runtime UI state in avivatorish stores.
- **Root Viv Config** — MDV-only legacy path where `VivMdvReact` serializes `config.viv.channelsStore` on the chart root. Not used by the SpatialData chart; not an upstream SpatialData.js concern.
- **App Viv Extensions** — runtime attachment: host app supplies Viv `LayerExtension` instances (and passes related props into the renderer). Extension classes are not serialized; extension data may live on Layer Channel Config when persistence is needed.
