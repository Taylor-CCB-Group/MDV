# Playwright Prompt And Context Guidance

This note records the issues encountered while adding a Playwright test for creating a 2D Scatter Plot, saving the view, and verifying persistence after reload. It also proposes what should be included in future prompts and in repo context files so the same issues do not recur.

## Issues Encountered

### 1. The documented Playwright command path did not run cleanly

The repo guidance pointed to:

```bash
pnpm run playwright-test ...
```

In this workspace, that script resolves to:

```bash
node node_modules/playwright/cli.js test
```

The checked-in tests import from `@playwright/test`. Under the current local setup, invoking the `playwright` CLI led to runner errors such as:

- `Playwright Test did not expect test.describe.configure() to be called here`
- `Requiring @playwright/test second time`

The test ran successfully only when invoked through the local `@playwright/test` CLI directly:

```bash
node node_modules/@playwright/test/cli.js test tests_playwright/project/scatter_view_persistence.spec.ts --project=chromium --reporter=list
```

### 2. Sandboxed browser launch failed

Running Chromium inside the sandbox failed with macOS permission errors during browser startup. The test needed an unsandboxed run to launch the browser.

This means any prompt asking the agent to both write and run Playwright tests should expect a permission escalation step for browser execution.

### 3. Backend availability was not the first real blocker

The initial sandboxed connectivity check to `http://localhost:5055` failed, but the real runtime blocker was browser-launch permission, not the backend itself. Once run outside the sandbox, the backend was reachable.

This means prompts should state clearly whether the backend is already running, but the context should also warn that sandboxed network checks can be misleading.

### 4. The Add Chart dialog needed a title before the chart could be created

Selecting `2D Scatter Plot` alone was not enough. The dialog left `Add Chart` disabled until a title was filled in.

The final test needed to:

- open the Add Chart dialog
- select `2D Scatter Plot`
- fill a unique chart title
- click `Add Chart`

### 5. The most stable persistence assertion was config-based, not purely visual

A brittle first attempt assumed the chart type would always be `wgl_scatter_plot`. That was too specific for this setup. The stable assertion was:

- create the chart with a unique title
- save the view
- reload
- verify the same chart config still exists after reload by `id`, `title`, and `type`

This is a better pattern for project-view persistence tests than relying only on chart counts or visuals.

### 6. Some existing project Playwright suites depend on stale setup helpers

Running the existing suites:

```bash
node node_modules/@playwright/test/cli.js test tests_playwright/project/chart_creation.spec.ts tests_playwright/project/view_management.spec.ts --project=chromium --reporter=list
```

failed before any test body ran.

The shared setup in `tests_playwright/utils/testUtils.ts`:

- calls `mockApiRoot(page)`
- navigates to `/`
- expects a visible `Create new project` control

In the current app state, the catalog page showed a populated `Recent Projects` view, and that exact `Create new project` text was not present. Both suites failed at:

- `tests_playwright/utils/testUtils.ts:13`

This is not a browser-launch issue and not a failure in the new scatter persistence test. It is a stale-assumption problem in the older project-test setup path.

### 7. There are now two distinct project-test patterns in the repo

The repo currently mixes:

- an older `newProjectSetup()` path in `tests_playwright/utils/testUtils.ts`
- a newer `createTemporaryProject()` path in `tests_playwright/utils/projectFixtures`

The newer helper is materially more reliable for backend-backed project tests because it:

- creates isolated temporary projects
- supports synthetic generated data
- cleans up after itself
- works against the real project page directly

The older helper depends on catalog-page assumptions and is now the source of immediate failures in the existing `chart_creation` and `view_management` suites.

## What Future Prompts Should Provide

Future prompts should explicitly provide the following information:

### Required prompt details

- The exact user flow to automate.
- Whether the backend is already running, and at which URL.
- Whether the agent is expected to run the tests after writing them.
- Whether the agent may request permission escalation for browser launch.
- Which helper or fixture path must be used for project creation.
- Whether older helpers such as `tests_playwright/utils/testUtils.ts` should be avoided.
- Which existing tests should be used as style references.
- What success condition matters most: UI presence, saved config persistence, reload persistence, or all three.

### Strongly recommended prompt details

- The exact Playwright command to use if the repo has a known CLI quirk.
- Whether the test should target a specific browser project such as `chromium`.
- Any UI prerequisites that are easy to miss, such as required chart titles or view names.
- Whether cleanup must happen through existing helpers.

## What Context Or Agent Reference Files Should Say Explicitly

The repo context files should make the following points explicit.

### In `AGENTS.md`

- The preferred Playwright execution command.
- The fallback command if `pnpm run playwright-test ...` fails because of runner/import mismatches.
- A note that browser launch in sandboxed agent sessions commonly needs escalation.
- A note that backend availability checks from inside the sandbox may not reflect the real host state.
- A note that backend-backed `project/` tests should prefer `tests_playwright/utils/projectFixtures` over `tests_playwright/utils/testUtils.ts` unless the task is explicitly about catalog flow.

Suggested addition:

```md
- for Playwright execution, prefer `pnpm run playwright-test ...`, but if the local runner errors with `@playwright/test` / `playwright` CLI mismatch messages, fall back to `node node_modules/@playwright/test/cli.js test ...`
- in sandboxed agent sessions, Chromium/WebKit/Firefox launch may require approval even when the backend is already running
- sandboxed localhost connectivity checks can be misleading; prefer confirming backend state from an unsandboxed run if browser launch is already being escalated
- for backend-backed tests under `tests_playwright/project/`, prefer `tests_playwright/utils/projectFixtures`; `tests_playwright/utils/testUtils.ts` follows an older catalog/setup path and may not match the current UI
```

### In `docs/PLAYWRIGHT_AGENT_WORKFLOW.md`

- Document the exact failure mode seen with the current script and import combination.
- Document the fallback command using `node node_modules/@playwright/test/cli.js test ...`.
- Document that some dialogs require additional inputs before primary actions become enabled.
- Document that for persistence tests, config-based assertions after reload are preferred over raw chart counts.
- Document that the older `testUtils.ts` helper is no longer the preferred setup path for backend-backed project tests.
- Document that if a test needs a real project page, the agent should start from `createTemporaryProject()` rather than from the catalog page.

Suggested addition:

```md
If `pnpm run playwright-test ...` fails with errors about `test.describe()` or duplicate `@playwright/test` loading, run:

node node_modules/@playwright/test/cli.js test <path-to-spec> --project=chromium --reporter=list

For Add Chart flows, do not assume selecting a chart type enables submission. Some chart types also require a title.

For view-persistence tests, prefer verifying the persisted chart config after reload using a unique chart title or id.

For backend-backed `project/` tests, prefer `tests_playwright/utils/projectFixtures`. Avoid `tests_playwright/utils/testUtils.ts` unless the test intentionally exercises the catalog-page project-creation UI.
```

### In `tests_playwright/README.md`

Add a short migration note so contributors do not copy the stale setup pattern into new tests.

Suggested addition:

```md
There are two setup styles in this folder. For backend-backed `project/` tests, prefer `utils/projectFixtures`.

Avoid using `utils/testUtils.ts` for new `project/` tests unless you specifically need to exercise the catalog page. That helper depends on older catalog UI assumptions and may fail before the project page is reached.
```

## Revised Prompt Template

Use the following prompt for similar work.

```md
Add a Playwright test for the project view flow that creates a 2D Scatter Plot, saves the current view, reloads the page, and verifies the chart persists.

Requirements:
- Use `tests_playwright/utils/projectFixtures` to create the project with synthetic data.
- Follow the existing Playwright patterns in `tests_playwright/project/`.
- Write the test in the `tests_playwright/project/` suite.
- Use a unique chart title so the persisted chart can be identified reliably after reload.
- Verify persistence using app/chart config after reload, not only visual presence.
- Clean up through the existing temporary-project helper.
- Do not use `tests_playwright/utils/testUtils.ts` for this flow.

Environment:
- The backend is already running at `http://localhost:5055`.
- The browser may need unsandboxed execution approval; request it if needed.
- First try the repo Playwright command from the repo root:
  `pnpm run playwright-test <spec> --project=chromium --reporter=list`
- If that fails with `@playwright/test` / `playwright` runner mismatch errors, fall back to:
  `node node_modules/@playwright/test/cli.js test <spec> --project=chromium --reporter=list`

Deliverables:
- The new Playwright spec.
- The command actually used to run it.
- The pass/fail result and any repo-level execution issue that should be documented.
```

## Shorter Final Prompt

If a shorter version is needed, use this:

```md
Write and run a Playwright test for creating a 2D Scatter Plot in a temporary synthetic-data project, saving the view, reloading, and verifying the saved chart persists.

Use `tests_playwright/utils/projectFixtures` and follow the existing `tests_playwright/project/` patterns. Give the chart a unique title and verify persistence from chart config after reload.

Do not use `tests_playwright/utils/testUtils.ts` for this flow.

The backend is already running at `http://localhost:5055`. Run from the repo root. Start with:
`pnpm run playwright-test <spec> --project=chromium --reporter=list`

If that fails with Playwright runner/import mismatch errors, use:
`node node_modules/@playwright/test/cli.js test <spec> --project=chromium --reporter=list`

If browser launch needs approval, request it and continue.
```

## Recommendation

For this repo, the most useful immediate improvement is to update the Playwright workflow documentation and `AGENTS.md` so they acknowledge the runner mismatch fallback explicitly. That will remove the biggest source of avoidable agent confusion on similar tasks.

The next most useful improvement is to mark `tests_playwright/utils/testUtils.ts` as a legacy setup path for project tests, or to migrate the older `project/` suites onto `createTemporaryProject()`.

## Existing Test Run Result

Using the working runner path:

```bash
node node_modules/@playwright/test/cli.js test tests_playwright/project/chart_creation.spec.ts tests_playwright/project/view_management.spec.ts --project=chromium --reporter=list
```

the result was:

- `2 failed`
- `25 did not run`

The failures were both in shared setup, before test execution:

- `tests_playwright/project/chart_creation.spec.ts`
- `tests_playwright/project/view_management.spec.ts`

Immediate cause:

- `tests_playwright/utils/testUtils.ts` expected `Create new project` on the catalog page
- the current page instead rendered the populated `Recent Projects` UI without that exact control

Recommended repo changes based on this run:

- update `AGENTS.md` with the fallback runner command and the preference for `projectFixtures`
- update `docs/PLAYWRIGHT_AGENT_WORKFLOW.md` with the same guidance
- update `tests_playwright/README.md` to mark `testUtils.ts` as a legacy project-test setup path
