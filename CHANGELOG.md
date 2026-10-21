# Changelog

All notable changes to `mdvtools` are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.4.0] - 2026-10-21

> ⚠️ Behaviour changes for the full application (`mdvtools[app]`): the server now uses
> SQLite by default and stores the database at `/app/mdv/mdv.sqlite3` unless you set
> `SQLITE_DB_PATH`. To keep using PostgreSQL, set `DB_BACKEND=postgres` before you start the
> server (#435). The server also records usage events, such as who opened which project and
> view, in a new `usage_events` table. Set `ENABLE_USAGE_TRACKING=0` to turn this off (#540).

A feature release. SQLite is now the default database for the full application, a new `SpatialData` chart overlays spatial layers, new charts show genomic data, and MDV works behind proxies that strip a path prefix, such as Verily Workbench. `pip install mdvtools` still installs the slim core.

### Added
- SQLite database backend for the full application, chosen with `DB_BACKEND` (`sqlite` or
  `postgres`). `deploy.sh`, `deploy.bat` and `deploy_gui.sh` now set up a local Docker deployment interactively (#435).
- A SpatialData chart, built on SpatialData.js, that overlays shapes, labels, points and image
  layers and has its own layers dialog (#514). Shapes, labels and points can be coloured and
  filtered from their associated table, and points layers have a feature list with
  per-feature selection and colour (#535, #551).
- `mdvtools convert-spatial --table-handling` merges multimodal tables, groups them by region
  or keeps them separate, and `--point-transform` sets how point coordinates are transformed
  (#530).
- Genome support: a `genome` datastore parameter, an SV Circos plot, a UCSC browser view, an
  embedded IGV browser and a category heat map. `mdvtools serve --track-dir` adds extra track
  folders. `pip install mdvtools` now installs `pysam` and `pyranges` for these (#470).
- A date column type. Datetimes from pandas and polars show as `YYYY-MM-DD` in tables,
  tooltips, axes, legends and filters, and still sort and filter in date order (#527).
- Legend filtering: click a category in a colour legend to filter by it, or drag across a
  continuous legend to filter a numeric range (#482, #493).
- ChatMDV can use a local model through Ollama or a choice of remote models, and you can pick
  which datasources a query uses (#513, #521).
- `mdvtools create-default-views` generates template views for a project (#548).
- `MDVProject.rename_column` and `MDVProject.soft_delete_column`, with matching project
  routes, so scripts can rename and hide columns (#547).
- Server extensions installed as separate Python packages, found through the
  `mdv.extensions` entry point and enabled by name in config (#540).
- Scripts in `mdvtools/scripts` for cleaning project data: `column_name_report.py` groups
  field names that differ only by case or separator, and `project_json_history.py` snapshots
  and restores `datasources.json`, `views.json` and `state.json` (#557).
- `VITE_SHOW_PROJECT_CREATION=true` at build time shows the new project and import project
  cards in a production build (#570).
- When a host frames MDV without cross-origin isolation, the catalog and project pages offer
  to open MDV in a new tab, where data loading works (#573).
- The debug dialog lists validation issues in a readable form (#491).

### Changed
- React Compiler is now enabled by default for Vite builds; set
  `VITE_USE_REACT_COMPILER=0` or `VITE_USE_REACT_COMPILER=false` to opt out locally (#502).
- The database assigns Project IDs, and new project folders get a random `uuid4` name instead
  of a name built from the ID. Projects copied into the projects folder get an owner, and
  rescans grant new projects to every admin (#543, #546).
- The server shows projects it cannot write as read-only in the interface. It used to serve
  them read-only without saying so (#544).
- Image layers support up to 10 channels after the viv update, and channel histograms are
  computed in a worker. MDV warns when an image's data layout will load slowly (#490, #565).
- The fraction legend is now a React component (#503).
- `sync_users_to_db` fetches the Auth0 admin role members once, instead of once per user
  (#526).
- Dependency updates: `@spatialdata/*` 0.12.0, deck.gl and luma.gl 9.4, viv 0.23 (#529, #541,
  #561).

### Fixed
- `pip install mdvtools` now installs `leidenalg` and `igraph`, so `--compute-x-umap` and
  `convert_scanpy_to_mdv(compute_x_umap=True)` can compute Leiden clusters. They were only in
  the `dev` dependency group before (#575).
- File upload, Socket.IO and page assets work behind a proxy that strips a path prefix
  (#571, #572).
- A user signed in to one deployment no longer carries their session into another
  deployment on the same host (#549).
- Spatial tables built by `by-region` or `merge` conversion colour and filter the right rows
  (#551).
- The points feature filter list stays responsive with thousands of features (#542).
- SpatialData shapes load in production builds (#539).
- Slow initial views no longer block the page behind a loading overlay (#519).
- Moving or resizing a chart no longer reorders the other charts, and dialogs can come to
  the front (#520).
- Box plot and violin plot x-axis labels update when their parameters change (#511), and dot
  plot x-axis labels fit the space available (#560).
- The gates column loads lazily, which avoids an out-of-memory crash on large datasets
  (#495).
- Popout tables work for both table types (#534).
- The user profile shows an error when it cannot load, instead of placeholder data (#489).
- Uploading to a new project no longer fails the permission check before its data file
  exists (#501).
- The server sends the headers that range requests need over CORS (#528).
- A rows-as-columns link can use a `unique` column as its name column (#533).

## [1.3.0] - 2026-06-04

> ⚠️ **Behaviour change — `mdvtools convert-spatial`** now takes the path to a SpatialData
> zarr store **directly**. It previously expected a *parent folder* containing one or more
> stores; pass the new `--batch` flag to restore that directory-scanning behaviour. Update
> any scripts or commands that pointed at a parent folder (#412).

A feature release: a new React data table, gating, new chart types, ChatMDV evaluation
tooling, and a large batch of fixes and dependency upgrades. Packaging was consolidated
internally — the published `pip install mdvtools` is the same slim/lite core as before.

### Added
- New **SlickGrid-based React data table**: inline cell editing, find-and-replace,
  multi-row highlighting, and sortable columns (#328).
- **Gating** — draw rectangle/polygon gates on scatter and spatial (viv) plots to define,
  name, colour, and manage cell populations (gates become a filterable column), with a
  Manage Gates dialog (#328, #388).
- Table column management: add, bulk-edit (fill all / fill empty), remove (impact dialog
  and soft delete), rename, and **compound columns — build a new column by combining two or
  more existing text columns** (#386, #436, #478).
- New chart types: a deck-based scatter **density** chart with shared contour settings and
  log bandwidth control (#390), and an experimental **Splatter Plot** (#391).
- **Permission shown in the GUI** — a view/edit lock icon on the project view (#359).
- Conversions: X-based UMAP & Leiden options for AnnData/SpatialData (#392); a helper to
  merge spatial (Xenium) annotations, patching `cell_id` (#444); a SpatialData conversion
  report runner (#412); tsv/tab file upload with file-extension hints (#405, #408).
- App/UX: Settings dialog with a General folder and search (#382); automatic view creation
  when a datasource is uploaded (#362); Enter-key submission in view dialogs (#377); React
  Text Box with collapsible markdown and mermaid diagrams (#425); subgroup selection in the
  link UI (#457); background filters (#460); schema-validation logging UI (#384); and
  category-selection widgets plus improved multitext/tag-annotation handling (#399, #401).
- ChatMDV: a verification step (#419) and an evaluation & testing framework (#469).
- The full application is installable from PyPI as an opt-in extra:
  `pip install "mdvtools[app]"` (database/server + chat/LLM + auth) (#477).

### Changed
- Consolidated the previously separate `pyproject.toml` files (the lite build under
  `mdvtools_lite_build/` and the full/dev project) into a single `python/pyproject.toml`.
  `pip install mdvtools` continues to install the same slim/lite core; the heavy
  server/chat/auth stack now lives behind the `app` extra, guarded at import (#477).
- Build and publish migrated from Poetry/Hatchling + twine to **uv** (#472).
- Supported Python is now `>=3.11,<3.13`.
- **Color Scheme dialog reimplemented in React and renamed "Color Palette"** (#415); the
  color legend likewise moved to a React component, deprecating the misspelled
  `overideValues` config key in favour of the correct spelling (backward compatible) (#479).
- Frontend toolchain: upgraded to Vite 8 (#442), migrated to pnpm (#441), CI to Node 24
  (#455); added **opt-in** React Compiler support to the Vite build (off by default; enable
  with `VITE_USE_React_COMPILER=1`) (#443). Dependency upgrades (#458); Playwright workflow
  tooling and synthetic SpatialData test data (#449, #432, #464, #473).
- ChartManager now boots through a React wrapper with a nicer loading state and more lazy
  module loading for bundle splitting (#448).
- Home navigation now uses links and preserves explicit directory routing (#389).

### Fixed
- Scatter rendering: local filter ownership, common grey layer, z-index/mouse events, and
  hiding missing values on grey layers (#427, #428, #429, #402).
- slickgrid row selection and stray React roots (#445); deck mouse-event rebinding for new
  versions (#447); React-chart class layout/overflow (#431); chart loading checks and
  clearer error messages (#465).
- Safari TLS error on the project view (#439); local dev and netlify preview routing
  (#385); invisible MUI modals in datasource fullscreen (#396).
- `add_datasource_polars` `TablePlot` now respects `supplied_columns_only` (#467);
  add-chart param state for multi-column selections (#400); React table replace now
  replaces the actual column value (#463); settings folder collapse during filtered
  search (#404).
- Wider scatter brush handles so they no longer disappear and are easier to click (#411).

[Unreleased]: https://github.com/Taylor-CCB-Group/MDV/compare/v1.4.0...HEAD
[1.4.0]: https://github.com/Taylor-CCB-Group/MDV/compare/v1.3.0...v1.4.0
[1.3.0]: https://github.com/Taylor-CCB-Group/MDV/releases/tag/v1.3.0
