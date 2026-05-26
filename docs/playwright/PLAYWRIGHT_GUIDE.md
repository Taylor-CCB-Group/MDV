# Playwright Guide

Developer reference for setting up, running, writing, and debugging Playwright tests.

**LLM / agent?** Read `docs/playwright/PLAYWRIGHT_AGENT.md` instead.

Reference docs:
- `PLAYWRIGHT_ARCHITECTURE.md` — fixtures, project lifecycle, execution policy
- `PLAYWRIGHT_WRITING_TESTS.md` — how to add and extend tests
- `PLAYWRIGHT_HISTORY.md` — background and prior failure modes
- `PLAYWRIGHT_ROADMAP.md` — planned work (CI, coverage, speed, playwright-cli)

---

## Quick start

Three things must be true before any Playwright work.

### 1. Docker backend

```bash
docker compose -f docker-secrets.yml up -d
curl -f http://localhost:5055
```

If the image is stale: `docker compose -f docker-secrets.yml up -d --build`

### 2. Browser binaries

```bash
pnpm run playwright-setup
```

Runs `pnpm install` and `pnpm exec playwright install --with-deps`. Run once after cloning
and again after any `pnpm install` that may have updated Playwright.

> **Do not** set `PLAYWRIGHT_BROWSERS_PATH` or install browsers to a custom path.
> Use Playwright's default browser location only.

### 3. Python environment (project tests only)

```bash
pnpm run playwright-preflight-project
```

Verifies Python interpreter, `h5py`/`scanpy` imports, `~/mdv` directory, and backend
reachability. If it fails, run `pnpm run python-setup` first.

---

## Running tests

Default target is `http://localhost:5055` (Docker stack). All project tests need the
backend running first.

### Project tests (backend-backed)

```bash
TEST_BASE_URL=http://localhost:5055 pnpm run playwright-test-project
```

Single spec:

```bash
TEST_BASE_URL=http://localhost:5055 pnpm run playwright-test-project tests_playwright/project/charts/chart_creation.spec.ts --project=chromium --reporter=list
```

### Catalog tests (browser-only, no backend lifecycle)

```bash
TEST_BASE_URL=http://localhost:5055 pnpm run playwright-test-catalog
```

### All tests

```bash
TEST_BASE_URL=http://localhost:5055 pnpm run playwright-test-all
```

### Playwright UI mode

```bash
TEST_BASE_URL=http://localhost:5055 pnpm run playwright-test-project-ui
TEST_BASE_URL=http://localhost:5055 pnpm run playwright-test-catalog-ui
TEST_BASE_URL=http://localhost:5055 pnpm run playwright-test-ui
```

Single spec in UI mode: append `--ui` to the base runner.

### Worktree / dev branch UI (5170 only)

Use only when testing UI changes not yet in the Docker image:

```bash
pnpm run dev    # separate terminal — Vite proxies API to localhost:5055
TEST_BASE_URL=http://127.0.0.1:5170 pnpm run playwright-test-project <spec> --project=chromium --reporter=list
```

### Suppress the HTML report

```bash
PLAYWRIGHT_HTML_OPEN=never TEST_BASE_URL=http://localhost:5055 pnpm run playwright-test-project ...
```

---

## Writing a new test

1. Copy `tests_playwright/project/_template.spec.ts` to `tests_playwright/project/<area>/<name>.spec.ts`.
2. Use `createTemporaryProjectViaSyntheticAnndata` — do not use `testUtils.ts` (legacy).
3. Wrap the test body in `try / finally { await handle.cleanup() }`.
4. Use helpers from `tests_playwright/utils/helpers.ts` for common UI actions.
5. Run with `pnpm run playwright-test-project <file> --project=chromium --reporter=list`.

For multiple related read-mostly tests sharing one project, see `chart_creation.spec.ts`
(`createSharedSyntheticAnndataSuite` + `beforeAll` / `afterAll`).

See `PLAYWRIGHT_WRITING_TESTS.md` for full guidance on catalog vs project tests, lifecycle
models, and anti-patterns.

---

## Exploring the UI with `playwright-cli`

`playwright-cli` is a browser exploration tool — not a test runner.

```bash
# Open a project page
pnpm exec playwright-cli open http://localhost:5055/project/<id> --headed

# Explore
pnpm exec playwright-cli snapshot
pnpm exec playwright-cli click e12
pnpm exec playwright-cli --raw generate-locator e12

# Optionally run a scratch script
pnpm exec playwright-cli run-code --filename=artifacts/playwright-cli/my-script.mjs

pnpm exec playwright-cli close
```

After exploring, write a fresh spec using the template. Do not promote `run-code` scripts
directly into tests — different syntax, no fixture system.

`artifacts/playwright-cli/` is gitignored scratch space.

---

## Debugging a failing test

Run with `--debug=cli` to pause the test and get a session name:

```bash
PLAYWRIGHT_HTML_OPEN=never pnpm run playwright-test-project tests_playwright/project/charts/chart_creation.spec.ts --project=chromium --reporter=list -- --debug=cli
```

Wait for `Debugging Instructions` and a session name (e.g. `tw-abcdef`), then attach in a
second terminal:

```bash
pnpm exec playwright-cli attach tw-abcdef
```

Take snapshots, inspect elements, step through the test. Every action emits equivalent
Playwright TypeScript you can copy into the spec.

---

## Cleanup

Fixture-created projects are cleaned up automatically. For projects from manual exploration:

```bash
pnpm run delete_projects            # dry-run by default; add --confirm to execute
PLAYWRIGHT_KEEP_PROJECTS=1          # env var to keep fixture projects while debugging
```

---

## Troubleshooting

**Missing browser binaries**

Run `pnpm run playwright-setup` from a normal terminal (not sandboxed). Do not install
browsers to a custom path — use Playwright's default browser location.

**Backend unreachable**

```bash
curl -f http://localhost:5055
docker compose -f docker-secrets.yml up -d
```

**Python env missing or broken**

```bash
pnpm run python-setup
pnpm run playwright-preflight-project --diagnostic
```

**Poetry / `.venv` blocked by macOS ACL metadata**

```bash
ls -le python/.venv
chmod -N python/.venv
cd python && poetry install --with dev
```

**Reset local Docker database**

```bash
docker compose -f docker-secrets.yml down -v
docker compose -f docker-secrets.yml up -d
```
