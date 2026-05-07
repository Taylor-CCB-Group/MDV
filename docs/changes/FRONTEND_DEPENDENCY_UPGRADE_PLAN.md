# Frontend Dependency Upgrade Plan

This document tracks dependency upgrade work that may be split across multiple PRs.

## Current goals

- Reduce `pnpm audit` findings without destabilizing the deck.gl/viv/luma stack.
- Prepare and execute React 19 migration while keeping chart behavior stable.
- Keep `@visx/*` available for React-side chart work and future refactors.
- Reserve deck/viv/luma major alignment and extension refactors for a dedicated PR when the pending Viv release is available.

## Audit-first track (completed)

### Applied in current audit-first updates

- Upgraded direct dependencies:
  - `socket.io-client` -> `4.8.3`
  - `uuid` -> `11.1.1`
- Upgraded dev dependencies:
  - `vite` -> `8.0.11`
  - `@vitejs/plugin-react` -> `^6.0.1` (required for Vite 8 compatibility)
  - `postcss` -> `8.5.14`
  - `tailwindcss` -> `3.4.19`
  - `jsdoc` -> `4.0.5`
- Added `pnpm.overrides` for transitive audit fixes where practical:
  - `fast-xml-parser` -> `4.5.5`
  - `lodash` -> `4.18.1`
  - `socket.io-parser` -> `4.2.6`
  - `preact` -> `10.29.1`
  - `mdast-util-to-hast` -> `13.2.1`
  - `protocol-buffers-schema` -> `3.6.1`
  - `postcss` -> `8.5.14`
  - `yaml` -> `2.8.4`
- Notes:
  - Keep `react-syntax-highlighter` on `16.x` for newer transitive dependencies/security posture; continue using `@types/react-syntax-highlighter` even though it's out of date until upstream ships bundled types.

### Known remaining risk after this phase

- Only one `pnpm audit` finding remains (`moderate`), for `fast-xml-parser` `<5.7.0`.
- The transitive path is `@deck.gl-community/editable-layers` -> `@deck.gl/geo-layers` -> `@loaders.gl/xml` -> `fast-xml-parser`.
- We upgraded to `fast-xml-parser` `4.5.5` to remove higher-severity findings, but are intentionally not forcing `5.x` in this pass due to possible behavior/API differences against the current deck/loaders dependency chain.
- The final `fast-xml-parser` remediation is planned for the dedicated deck/viv/luma alignment PR below.

## Planned PR breakdown

### PR A: Audit-first non-breaking dependency updates

- Scope:
  - direct + transitive updates that do not require deck/viv/luma major coordination
  - lockfile refresh
  - confirm build and tests

### PR B: React 19 migration with visx retained

- Scope:
  - `react`/`react-dom` + type packages to `19.x`
  - update to latest `@visx/*` packages and verify compatibility on React 19
  - migrate from the current PR build of `slickgrid-react` to a cleaner current release compatible with the React 19 stack
  - update any React 19 type/ref adjustments
  - document known chart-specific regressions if discovered
- Notes:
  - Prior attempt removed visx and broke axis rendering; this should be explicitly avoided.

### PR C: Viv/deck/luma alignment after new Viv release

- Scope:
  - update `@hms-dbmi/viv` to the new published version
  - align `deck.gl`, `@deck.gl/*`, `@deck.gl-community/*`, `@luma.gl/*`, and `@loaders.gl/*`
  - remove temporary overrides that are no longer needed
  - resolve peer-dependency mismatches (especially community layers)
- Expected code work:
  - refactor extensions for UBO changes
  - adjust variable Viv channels integration
  - validate no duplicate luma versions at runtime

### PR D: React performance track (optional, follow-up)

- Scope:
  - profile React rendering hot paths in chart-heavy screens
  - evaluate React Compiler compatibility and rollout strategy
  - capture win/loss metrics before broad enablement

### PR E: Vite 8 readiness and rollout

- Context:
  - The earlier `RefreshRuntime is not defined` issue was caused by a version mismatch (`vite@8` with `@vitejs/plugin-react@5.x`, whose peer range only supports Vite 4-7).
  - Updating to `@vitejs/plugin-react@6.x` resolves the Vite 8 compatibility gap.
  - There is also prior memory of `vitest` compatibility issues around decorator-heavy codepaths (for example `src/links/link_utils.ts`), which should be treated as a first-class migration risk.
- Goal:
  - Keep Vite 8 if dev and test behavior remain stable, and capture potential build-time/tooling improvements.
- Scope:
  - validate `@vitejs/plugin-react` + `vite` + `vitest` version matrix for Vite 8
  - reproduce and fix decorator-related test/runtime incompatibilities
  - re-run dev server + HMR smoke tests for dialog/chart-heavy paths
  - compare build times (`pnpm run vite-build`) between the last Vite 7 baseline and Vite 8.x
- Exit criteria:
  - no React refresh runtime errors in dev
  - `vitest` suite passes without decorator regressions
  - no new `pnpm audit` findings introduced
  - measured build-time improvement or a clear justification to defer

## Validation checklist for each dependency PR

- `pnpm install --no-frozen-lockfile`
- `pnpm audit`
- `pnpm run test`
- `pnpm run vite-build`
- targeted chart smoke checks (deck + viv + visx-backed components)
