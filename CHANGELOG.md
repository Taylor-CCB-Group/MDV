# Changelog

All notable changes to `mdvtools` are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.3.0] - 2026-06-04

A feature release: new chart types, table editing, gating/permissions, ChatMDV
evaluation tooling, and a large batch of fixes and dependency upgrades. Packaging was
consolidated internally — the published `pip install mdvtools` is the same slim/lite
core as before.

### Added
- Gating feature and access controls, plus permission display in the GUI (#328, #388, #389, #359).
- New chart types: DeckSplatter (#391) and deck scatter density with shared contour
  settings and log bandwidth control (#390).
- Table editing in the new react table: add, bulk-edit, remove, and rename columns, plus
  hard compound columns (#386, #436, #478).
- Conversions: X-based UMAP & Leiden options for AnnData/SpatialData (#392), Xenium
  annotation merge (#444), a SpatialData conversion report runner with `--batch` flag
  (#412), and tsv/tab file upload support (#405, #408).
- App/UX: Settings dialog with general folder and search (#382), view creation for a new
  datasource (#362), keyboard shortcuts for view dialogs (#377), react Text Box with
  collapsible markdown and mermaid diagrams (#425), subgroup selection in the link UI
  (#457), background filters (#460), react color scheme dialog (#415), schema-validation
  logging UI (#384), and category-selection widgets (#399, #401).
- ChatMDV: verification (#419) and an evaluation & testing framework (#469).
- The full application is now installable from PyPI as an opt-in extra:
  `pip install "mdvtools[app]"` (database/server + chat/LLM + auth).

### Changed
- Consolidated the previously separate `pyproject.toml` files (the lite build under
  `mdvtools_lite_build/` and the full/dev project) into a single `python/pyproject.toml`.
  `pip install mdvtools` continues to install the same slim/lite core; the heavy
  server/chat/auth stack now lives behind the `app` extra, guarded at import (#477).
- Build and publish migrated from Poetry/Hatchling + twine to **uv** (#472).
- Supported Python is now `>=3.11,<3.13`.
- Frontend toolchain: upgraded to Vite 8 (#442), migrated to pnpm (#441), enabled the
  React compiler (#443), and moved CI to Node 24 (#455); various dependency upgrades
  (#458). Playwright workflows stabilized with added tooling (#449, #432, #464, #473).
- ChartManager now boots through a react wrapper with a nicer loading state and more lazy
  module loading for bundle splitting (#448).

### Fixed
- Scatter rendering: local filter ownership, common grey layer, z-index/mouse events, and
  hiding missing values on grey layers (#427, #428, #429, #402).
- slickgrid row selection and stray react roots (#445); deck mouse-event rebinding for new
  versions (#447); react-chart class layout/overflow (#431); chart loading checks and
  clearer error messages (#465).
- Safari TLS error on the project view (#439); local dev and netlify preview routing
  (#385); invisible MUI modals in datasource fullscreen (#396).
- `add_datasource_polars` `TablePlot` now respects `supplied_columns_only` (#467);
  add-chart param state for multi-column selections (#400); react table replace now
  replaces the actual column value (#463); settings folder collapse during filtered
  search (#404).

[Unreleased]: https://github.com/Taylor-CCB-Group/MDV/compare/v1.3.0...HEAD
[1.3.0]: https://github.com/Taylor-CCB-Group/MDV/releases/tag/v1.3.0
