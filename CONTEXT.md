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

## Tooling terms (kept distinct — the proposal conflated them)

- **Build backend** — the tool that turns the source tree into a wheel/sdist
  (e.g. `uv_build`, `hatchling`, `poetry-core`). Determines how files are selected into
  the wheel.

- **Installer / lockfile manager** — the tool that resolves, locks, and installs
  dependencies into an environment (e.g. `uv`, `pip`, `poetry`). Produces the lockfile
  (`uv.lock`). Independent of the build backend.

## SpatialData.js integration terms

Canonical rendering vocabulary lives in [SpatialData.js CONTEXT.md](https://github.com/Taylor-CCB-Group/SpatialData.js/blob/main/CONTEXT.md). MDV chart code uses the same terms:

- **Render Stack** — saved draw order on chart config (`renderStack.entries`), not parallel `layerOrder` arrays
- **Stack Entry** — one ordered item (`kind: spatial | host | group`)
- **Host Overlay** — MDV deck layer referenced by `hostLayerId` in the stack, resolved at runtime
- **Runtime Attachment** — `hostLayerResolver`, tooltip hooks, `deckProps`; not serialized in `entry.props`
- **MobX Control Island** — layer dialog UI that patches `config.renderStack` directly

See [docs/spatialdata-vis-integration.md](docs/spatialdata-vis-integration.md) for MDV-specific integration notes.