# Changelog

All notable changes to `mdvtools` are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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

[Unreleased]: https://github.com/Taylor-CCB-Group/MDV/compare/v1.3.0...HEAD
[1.3.0]: https://github.com/Taylor-CCB-Group/MDV/releases/tag/v1.3.0
