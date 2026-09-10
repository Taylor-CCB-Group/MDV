# MDV Context Glossary

Canonical terms for MDV. Glossary only: no implementation detail, no decisions.
Decisions that are hard to reverse live in `docs/adr/`.

## Packaging terms

- **mdvtools**: the single Python package published to PyPI (`pip install mdvtools`).
  There is one package and one version; "slim" and "full" are *consumption modes* of it,
  not separate products.

- **Slim core**: the **required** dependency set: what a bare `pip install mdvtools`
  gives you. Targeted at the majority of users (manipulating / viewing MDV projects).
  Includes the spatial feature (spatialdata is a required dependency, not an extra).

- **Extra**: an optional dependency group declared under
  `[project.optional-dependencies]`, installed on request via `mdvtools[<name>]`. Heavy /
  minority features only.

- **`app` extra**: the single extra. Pulls the entire optional cluster: database/server
  (sqlalchemy, psycopg2, gunicorn), chat/LLM (langchain, faiss, matplotlib), and auth
  (authlib, jose, auth0, redis). These three are a runtime import *cycle*
  (`auth → dbutils → llm → auth`), so they cannot be installed independently; they are
  one bundle. `mdvtools[app]` == the full app; there is no coherent middle ground between
  slim core and `app`.

- **Full install**: slim core **plus the `app` extra**. What the root `Dockerfile` builds
  and runs (database, view gallery, chat, auth). Never the default for a PyPI user.

- **Guarded import**: a top-level heavy import in a leaf module wrapped so that, on a slim
  install, calling the feature raises a friendly "install `mdvtools[<extra>]`" error
  instead of an `ImportError` at import time.

## Tooling terms (kept distinct: the proposal conflated them)

- **Build backend**: the tool that turns the source tree into a wheel/sdist
  (e.g. `uv_build`, `hatchling`, `poetry-core`). Determines how files are selected into
  the wheel.

- **Installer / lockfile manager**: the tool that resolves, locks, and installs
  dependencies into an environment (e.g. `uv`, `pip`, `poetry`). Produces the lockfile
  (`uv.lock`). Independent of the build backend.

## Database-backed catalog terms

**Independent MDV deployment**:
One MDV application with its own PostgreSQL or SQLite database, lifecycle, Project IDs,
permissions, and endpoint. Two of them may share a Project root without sharing anything else.
_Avoid_: replica, copy, instance

**Project ID**:
The database-local integer identity of a project in one deployment's catalog.
_Avoid_: directory name, folder ID, UUID, global project identity

**Project root**:
The configured filesystem root that contains project directories. May be shared between
Independent MDV deployments, by mount or by bucket sync.
_Avoid_: project ID, catalog

**Recovery copy**:
The fields a rescan needs to rebuild a catalog row, written into the project directory
(`state.json`): the display name and the permission. Read only when no row exists.
_Avoid_: sidecar, source of truth, project metadata

**Route registry**:
A process-local, disposable projection of active projects used to dispatch project requests.
Rebuilt from the catalog at startup.
_Avoid_: catalog, source of truth, shared routing table

**Shared project**:
Project files visible to more than one Independent MDV deployment.
_Avoid_: replica storage, collaborative project, globally identified project

**Storage name**:
The directory name a project occupies under the Project root. Never parsed and never derived
from the Project ID. New projects get a random one; directories that already exist keep the
names they have.
_Avoid_: project folder ID, project slug, directory ID

**Unowned project**:
A project with a catalog row and no `user_projects` row marking an owner. A project copied
into the Project root arrives this way, because the scan that finds it has no user to
attribute it to. Startup and rescan give every administrator ownership of one, since the
project list only shows a user the projects they hold a permission row for.
_Avoid_: orphan, public project, shared project

## SpatialData.js integration terms

Canonical rendering vocabulary lives in [SpatialData.js CONTEXT.md](https://github.com/Taylor-CCB-Group/SpatialData.js/blob/main/CONTEXT.md). MDV chart code uses the same terms:

- **Render Stack**: saved draw order on chart config (`renderStack.entries`), not parallel `layerOrder` arrays
- **Stack Entry**: one ordered item (`kind: spatial | host | group`)
- **Host Overlay**: MDV deck layer referenced by `hostLayerId` in the stack, resolved at runtime
- **Runtime Attachment**: `hostLayerResolver`, tooltip hooks, `deckProps`; not serialized in `entry.props`
- **Image Layer Registry**: chart-owned runtime attachment (`SpatialDataMdvReact`): callbacks that expose loaded Viv image data (`getImageLoadedDataByElementKey`, load state) from the spatial renderer. The viewer populates it; the layer dialog consumes it via `ImageLayerContextProvider`. Bridges separate React trees (chart vs dialog portal) without duplicating image load paths.
- **Spatial Image Panel Context**: React context on [`ImageLayerPanel`](src/react/components/spatialLayers/ImageLayerPanel.tsx): wraps `ImageLayerContextProvider` (upstream loaded-image defaults) and exposes `useLayerChannelState` write API (`setChannels`, `addChannel`, `removeChannel`) to channel UI. Persisted field edits go through the hook, not zustand `channelsStore`.
- **Image Layer Runtime Bridge**: pure helpers in `image_layer_runtime.ts`: `channelId`-keyed stats cache (`domains`, `raster`), viewer parallel-array sync (`channelOptions`, loading flags), and tone to `vivLayerProps` helpers. Not a second MobX persistence bridge.
- **MobX Control Island**: layer dialog UI that patches `config.renderStack` directly
- **Layer Channel Config**: serializable image channel state on `renderStack.entries[].props.channels` (`ChannelConfig` in `@spatialdata/vis`): colors, contrast limits, visibility, selections, and extension-related data fields. Tone (brightness/contrast) is **not** part of `channels`; it persists in `entry.props.vivLayerProps`. Distinct from runtime UI state in avivatorish stores.
- **Root Viv Config**: MDV-only legacy path where `VivMdvReact` serializes `config.viv.channelsStore` on the chart root. Not used by the SpatialData chart; not an upstream SpatialData.js concern.
- **App Viv Extensions**: runtime attachment: host app supplies Viv `LayerExtension` instances (and passes related props into the renderer). Extension classes are not serialized; extension data may live on Layer Channel Config when persistence is needed.
