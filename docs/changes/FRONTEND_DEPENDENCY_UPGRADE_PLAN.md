# Frontend Dependency Upgrade Plan

This document tracks dependency upgrade work that may be split across multiple PRs.

## Current goals

- Reduce `pnpm audit` findings without destabilizing the deck.gl/viv/luma stack.
- Keep the React 19 migration stable while chart behavior is verified.
- Keep `@visx/*` available for React-side chart work and future refactors.
- Reserve deck/viv/luma major alignment, SpatialData.js integration, and extension refactors for dedicated follow-up work.

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
  - `lodash` -> `4.18.1`
  - `socket.io-parser` -> `4.2.6`
  - `preact` -> `10.29.1`
  - `mdast-util-to-hast` -> `13.2.1`
  - `protocol-buffers-schema` -> `3.6.1`
  - `postcss` -> `8.5.14`
  - `yaml` -> `2.8.4`
  - `brace-expansion` -> `5.0.6`
  - `fast-xml-parser` -> `5.7.0`
  - `ws` -> `8.20.1`
- Notes:
  - Keep `react-syntax-highlighter` on `16.x` for newer transitive dependencies/security posture; continue using `@types/react-syntax-highlighter` even though it's out of date until upstream ships bundled types.

### Known remaining risk after this phase

- The previous remaining `fast-xml-parser` audit item was remediated by overriding to `5.7.0`.
- Keep watching the deck/loaders path (`@deck.gl-community/editable-layers` -> `@deck.gl/geo-layers` -> `@loaders.gl/xml` -> `fast-xml-parser`) during the dedicated deck/viv/luma alignment PR because that stack may eventually remove the need for the override.

## Planned PR breakdown

### PR A: Audit-first non-breaking dependency updates (completed)

- Scope:
  - direct + transitive updates that do not require deck/viv/luma major coordination
  - lockfile refresh
  - confirm build and tests

### PR B: React 19 migration with visx retained (completed)

- Scope:
  - `react`/`react-dom` + type packages to `19.x`
  - update to latest `@visx/*` packages and verify compatibility on React 19
  - migrate from the current PR build of `slickgrid-react` to a cleaner current release compatible with the React 19 stack
  - update any React 19 type/ref adjustments
  - document known chart-specific regressions if discovered
- Notes:
  - Prior attempt removed visx and broke axis rendering; this should be explicitly avoided.
  - Execution notes (May 2026):
    - React and React DOM were moved to `19.2.6`, with `@types/react`/`@types/react-dom` updated to current `19.x`.
    - `slickgrid-react` and `@slickgrid-universal/common` were migrated from PR snapshot URLs to stable `10.6.0` releases.
    - `@visx/*` was moved to `4.0.0-alpha.11`, keeping visx-backed axis rendering available on the React 19 stack.
    - Additional React 19 peer-range warnings are also present for `mobx-react-lite`, `@react-spring/*` (via visx), `use-sync-external-store`, and `@welldone-software/why-did-you-render`.

### PR C: React performance track (compiler enabled by default)

- Scope:
  - keep React Compiler compatibility and rollout notes in one focused tooling PR
  - profile React rendering hot paths in chart-heavy screens
  - capture win/loss metrics before broad enablement
- Current intent (June 2026):
  - Keep runtime/app changes outside this track; those were already cherry-picked and do not require a separate PR in this worktree.
  - React Compiler is enabled by default in Vite builds.
  - Use `VITE_USE_REACT_COMPILER=0` or `VITE_USE_REACT_COMPILER=false` to opt out for local comparison or emergency fallback.
- Initial React Compiler enablement notes:
  - `@vitejs/plugin-react` v6 exposes React Compiler through `reactCompilerPreset()` and `@rolldown/plugin-babel`.
  - The compiler preset is currently filtered to `.tsx` files. The default preset filter also caught vanilla chart modules with decorators, which made the compiler Babel pass fail while parsing decorated chart classes.
  - Direct Vite production build validation passes with the `.tsx` filter:
    `pnpm exec cross-env NODE_OPTIONS="--max-old-space-size=4096" build=dev_pt vite build --outDir ./vite-dist`
  - First runtime gotcha found: vanilla code (e.g. `ChartManager.js`, `static_index.ts`, `DataLoaderUtil.ts`) calling TSX wrapper components as plain functions (`MenuBarWrapper()`, `ErrorComponentReactWrapper({...})`, `ProjectStateHandlerWrapper({...})`) triggers `Invalid hook call` because the compiler injects `useMemoCache` into compiled functions. Fix is to render via `React.createElement(Comp, props)` (or JSX) so the function runs inside React's render pipeline. Other similar invocations should be audited as part of broader rollout.
- Local smoke-test direction:
  - Use the catalog Playwright tests against a Vite dev server for React shell and route-mocked catalog workflows.
  - Use the project Playwright tests against a live backend for chart-manager, datastore, vanilla-event, MobX, and Zustand interaction coverage.
  - `docs/playwright/PLAYWRIGHT_GUIDE.md` documents the current commands; add a named PR C smoke subset there if this track becomes a regular local validation step.
- Latest local validation note (8th May 2026):
  - Off/on probe runs against `scripts/playwright/playwright_compiler_probe.mjs` completed in both modes without hard runtime failures.
  - Timing and payload composition varied between runs (different chart/data composition observed), so current measurements are treated as smoke validation, not a strict performance benchmark.
- Latest build validation note (9th June 2026):
  - `pnpm exec tsgo` passes.
  - Direct Vite production build validation passes with React Compiler enabled by default:
    `pnpm exec cross-env NODE_OPTIONS="--max-old-space-size=4096" build=dev_pt vite build --outDir ./vite-dist`
- React Compiler lint follow-up:
  - Added a dedicated probe command:
    `pnpm run lint:react-compiler`
  - The probe uses `eslint-plugin-react-hooks` `recommended-latest` rules through `eslint.react-compiler.config.mjs`.
  - Treat this as diagnostic output, not a CI gate yet.
  - Current scope is `src/**/*.{ts,tsx}`, excluding `src/react/components/avivatorish/**` because legacy inline ESLint comments there add noise unrelated to compiler readiness.
  - First run reports 143 diagnostics across 52 files:
    - `react-hooks/exhaustive-deps`: 52 warnings
    - `react-hooks/set-state-in-effect`: 38 errors
    - `react-hooks/refs`: 21 errors
    - `react-hooks/immutability`: 19 errors
    - `react-hooks/rules-of-hooks`: 9 errors
    - `react-hooks/use-memo`: 3 errors
    - `react-hooks/preserve-manual-memoization`: 1 error
  - Suggested triage order:
    - fix real hook ordering and dependency issues first
    - review ref reads during render and config/prop mutation in React components
    - convert derived state-in-effect patterns where practical
    - keep imperative integration exceptions explicit for chart, deck.gl, SlickGrid, MobX, and popout/window plumbing


### PR D: SpatialData.js integration after Viv/deck/luma alignment

- Scope:
  - package alignment for `@hms-dbmi/viv`, `deck.gl`, `@deck.gl/*`, `@deck.gl-community/*`, `@luma.gl/*`, and `@loaders.gl/*` is done in this worktree
  - continue WIP integration of the first released SpatialData.js packages (`@spatialdata/core@0.1.0` and `@spatialdata/react@0.1.0`) where they replace bespoke SpatialData store access
  - remove temporary overrides that are no longer needed
  - resolve peer-dependency mismatches (especially community layers)
- Expected code work:
  - refactor extensions for UBO changes
  - move WIP SpatialData.js integration from prototype paths into the supported Viv/spatial image flow
  - adjust variable Viv channels integration
  - validate no duplicate luma versions at runtime



## Validation checklist for each dependency PR

- `pnpm install --no-frozen-lockfile`
- `pnpm audit`
- `pnpm run test`
- `pnpm run vite-build`
- targeted chart smoke checks (deck + viv + visx-backed components)
