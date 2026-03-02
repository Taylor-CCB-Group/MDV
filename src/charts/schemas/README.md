# MDV Chart & Datasource Schemas

This directory contains Zod schemas for MDV chart configurations and datasources.  
They serve two main purposes:

- **Runtime validation** in the browser (logging when configs don’t conform).
- **Generating JSON Schema artefacts** for external tooling (treated as build outputs).

## Key files

- `ChartConfigSchema.ts`
  - Defines:
    - `BaseConfigSchema`
    - `ChartColorConfigSchema`
    - Per‑chart schemas (scatter, bar, histogram, etc.)
    - `ChartConfigSchema` (union of all chart config types)
  - Registers each chart schema with `ChartConfigRegistry` so we can:
    - Look up a schema by chart `type`/`version`.
    - Build the `ChartConfigSchema` union from the registry.

- `ChartConfigRegistry.ts`
  - Simple in‑memory registry:
    - `registerChartConfigSchema(type, schema, { version })`
    - `getChartConfigSchema(type, version?)`
    - `getAllChartConfigSchemas()`

- `DataSourceSchema.ts`
  - Zod schemas for datasource configs:
    - Column metadata (`DataSourceColumnSchema`)
    - Column groups, links, images, regions, interactions, offsets, genome browser, etc.
  - Exposes:
    - `DataSourceSchema`
    - `DataSourcesArraySchema`
    - `validate*/safeValidate*` helpers.

- `scripts/generate-schemas.ts`
  - Used by `npm run build-schemas`.
  - Converts the Zod schemas to JSON Schema and writes:
    - `schemas/chart-config-schema.json`
    - `scripts/schemas/chart-config-schema.json`
    - `scripts/schemas/datasource-schema.json`
    - `scripts/schemas/datasources-array-schema.json`
  - These JSON files are **generated artefacts** and are git‑ignored.

## How schemas are used at runtime

- **Chart configs**
  - `initialiseChartConfig` (in `src/charts/chartConfigUtils.ts`) calls:
    - `getChartConfigSchema(config.type, config.version)` and falls back to `BaseConfigSchema`.
    - `safeParse` on the chosen schema.
    - On failure, logs a validation error (and records it for the Debug JSON dialog).
  - Charts still render on a best‑effort basis; validation is log‑only in this PR.

- **Datasource configs**
  - `DataStore` (in `src/datastore/DataStore.js`) calls:
    - `DataSourceSchema.safeParse(config)` in the constructor.
    - On failure, logs the issues and records them for inspection in the Debug JSON dialog.
  - Behaviour is again log‑only; we keep using the original config.

## Generating JSON Schema

To regenerate JSON Schema files from the Zod definitions:

```bash
npm run build-schemas
```

This:

- Runs `scripts/generate-schemas.ts`.
- Emits JSON Schema files under `schemas/` and `scripts/schemas/`.
- Leaves the TypeScript Zod definitions here as the single source of truth.

The generated JSON files are intended for:

- External validation tooling.
- Potential future Python / API integrations.

They are **not** committed to the repo and can be safely deleted and regenerated.

## Future directions (informal)

Some directions we may explore, building on these schemas:

- **Schema‑driven UI** for chart creation/settings.
- Gradually moving from generic `param` arrays to explicit, named properties.
- Richer column‑type constraints encoded in the schemas.