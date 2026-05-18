/**
 * Backend project fixtures for Playwright `project/` specs.
 *
 * | File | What it does |
 * |------|----------------|
 * | `importZip.ts` | Mock zip → `POST /import_project` (+ optional CSV fallback) |
 * | `syntheticAnndata.ts` | `generate_synthetic_anndata_project.py` + `/rescan_projects` |
 * | `syntheticSpatial.ts` | `generate_synthetic_spatial_project.py` + `/rescan_projects` |
 * | `core.ts` | Shared: wait ready, delete API, list projects, rescan |
 * | `pythonEnv.ts` | Poetry / `PYTHONPATH` for subprocess helpers |
 */

export * from "./core";
export * from "./importZip";
export * from "./pythonEnv";
export * from "./syntheticAnndata";
export * from "./syntheticSpatial";
