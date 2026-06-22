# Playwright — Agent Reference

Compact rules and prompt templates for LLM/agent sessions.

**Human developer?** Read `docs/playwright/PLAYWRIGHT_GUIDE.md` instead.

Planned improvements: `docs/playwright/PLAYWRIGHT_ROADMAP.md`

---

## Non-negotiable rules

- Only run backend project tests in an **unsandboxed** environment (need `localhost:5055` and a real browser).
- Do **not** set `PLAYWRIGHT_BROWSERS_PATH` or install browsers to a custom/sandbox path. If browsers are missing, tell the user to run `pnpm run playwright-setup` in their terminal.
- Do **not** use `pnpm run playwright-test -- ...` — injects a literal `--` into Playwright's args.
- Do **not** use `pnpm exec playwright test` — wrong config, may hit npm registry.
- Do **not** use `tests_playwright/utils/testUtils.ts` for new project specs — it is legacy.
- `pnpm exec playwright-cli` is for **exploration only** — not a test runner.
- After `playwright-cli` exploration, always write a **fresh spec** from `_template.spec.ts`. Do not convert `run-code` scripts directly into tests.
- `artifacts/playwright-cli/` is gitignored scratch space — nothing there is a committed test.

---

## Before doing any Playwright work

Run all three in order:

```bash
# 1. Backend
docker compose -f docker-secrets.yml up -d
curl -f http://localhost:5055

# 2. Browser binaries (run in user terminal only)
pnpm run playwright-setup

# 3. Python + preflight (project tests only)
pnpm run playwright-preflight-project
```

If preflight fails, the user needs to run `pnpm run python-setup` first.

---

## Running tests

| Task | Command |
|------|---------|
| Project tests (default) | `TEST_BASE_URL=http://localhost:5055 pnpm run playwright-test-project` |
| Single spec | `TEST_BASE_URL=http://localhost:5055 pnpm run playwright-test-project <file> --project=chromium --reporter=list` |
| Catalog tests | `TEST_BASE_URL=http://localhost:5055 pnpm run playwright-test-catalog` |
| All tests | `TEST_BASE_URL=http://localhost:5055 pnpm run playwright-test-all` |
| Branch UI (Vite) | `TEST_BASE_URL=http://127.0.0.1:5170 pnpm run playwright-test-project <file> --project=chromium --reporter=list` |
| Suppress HTML report | Prefix with `PLAYWRIGHT_HTML_OPEN=never` |

---

## Writing a test

1. Copy `tests_playwright/project/_template.spec.ts` to `tests_playwright/project/<area>/<name>.spec.ts`.
2. Use `createTemporaryProjectViaSyntheticAnndata` — not `testUtils.ts`.
3. Wrap test body in `try / finally { await handle.cleanup() }`.
4. Use `tests_playwright/utils/helpers.ts` for common UI actions.
5. Validate: `TEST_BASE_URL=http://localhost:5055 pnpm run playwright-test-project <file> --project=chromium --reporter=list`

---

## `playwright-cli` exploration workflow

```bash
pnpm run playwright-preflight-project
pnpm exec playwright-cli open http://localhost:5055/project/<id> --headed
# snapshot / click / generate-locator as needed
pnpm exec playwright-cli close
```

Then write a fresh spec — do not convert `run-code` scripts directly.

---

## Prompt templates

Copy one into a new agent chat and fill in the bracketed parts.

### Write a new project test (Docker / 5055)

```
Read docs/playwright/PLAYWRIGHT_AGENT.md first.

Task: write a Playwright project test for [describe the flow].

Environment:
- Backend on http://localhost:5055 — run pnpm run playwright-preflight-project first
- TEST_BASE_URL=http://localhost:5055
- Unsandboxed environment required
- Do not set PLAYWRIGHT_BROWSERS_PATH or install browsers in sandbox

Implementation:
- Copy tests_playwright/project/_template.spec.ts
- Use createTemporaryProjectViaSyntheticAnndata only (not testUtils.ts)
- try/finally cleanup on the project handle
- Use helpers from tests_playwright/utils/helpers.ts

Run:
TEST_BASE_URL=http://localhost:5055 pnpm run playwright-test-project tests_playwright/project/<area>/<name>.spec.ts --project=chromium --reporter=list
```

### Write a test for branch UI (Vite / 5170)

```
Read docs/playwright/PLAYWRIGHT_AGENT.md first.

Task: write a Playwright project test for [describe the UI flow on this branch].

Environment:
- Backend on http://localhost:5055, frontend via pnpm run dev (http://127.0.0.1:5170)
- TEST_BASE_URL=http://127.0.0.1:5170
- Unsandboxed environment required
- Do not set PLAYWRIGHT_BROWSERS_PATH or install browsers in sandbox

Implementation:
- Copy tests_playwright/project/_template.spec.ts
- Use createTemporaryProjectViaSyntheticAnndata only
- try/finally cleanup on the project handle

Run:
TEST_BASE_URL=http://127.0.0.1:5170 pnpm run playwright-test-project tests_playwright/project/<area>/<name>.spec.ts --project=chromium --reporter=list
```

### Explore with playwright-cli, then write a spec

```
Read docs/playwright/PLAYWRIGHT_AGENT.md first.

Task: explore [describe the flow] with playwright-cli, then write a committed Playwright spec.

Explore (artifacts/playwright-cli/ is scratch — do not commit):
- pnpm run playwright-preflight-project
- pnpm exec playwright-cli open http://localhost:5055/project/<id> --headed
- snapshot / click / generate-locator as needed
- pnpm exec playwright-cli close

Then write a fresh spec:
- Copy tests_playwright/project/_template.spec.ts
- Use projectFixtures + helpers.ts
- Validate: pnpm run playwright-test-project <spec> --project=chromium --reporter=list

Do not convert run-code scripts directly into tests.
```

### Run or fix an existing spec

```
Read docs/playwright/PLAYWRIGHT_AGENT.md first.

Task: [run / fix] tests_playwright/project/<path>/<name>.spec.ts

Environment:
- TEST_BASE_URL=http://localhost:5055 (or http://127.0.0.1:5170 for branch UI on Vite)
- Run pnpm run playwright-preflight-project first
- Unsandboxed environment required
- Do not set PLAYWRIGHT_BROWSERS_PATH or install browsers in sandbox

Run:
PLAYWRIGHT_HTML_OPEN=never TEST_BASE_URL=<url> pnpm run playwright-test-project tests_playwright/project/<path>/<name>.spec.ts --project=chromium --reporter=list
```
