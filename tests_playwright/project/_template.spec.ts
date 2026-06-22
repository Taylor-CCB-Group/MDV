/**
 * Template for backend-backed project tests.
 *
 * Quick start:
 *   1. Copy to tests_playwright/project/<area>/<name>.spec.ts
 *   2. Replace "Template" with a descriptive suite name.
 *   3. Implement the test body.
 *   4. Run:
 *        TEST_BASE_URL=http://localhost:5055 pnpm run playwright-test-project tests_playwright/project/<area>/<name>.spec.ts --project=chromium --reporter=list
 *
 * Prerequisites (see docs/playwright/PLAYWRIGHT_GUIDE.md):
 *   - Docker backend running on localhost:5055
 *   - `pnpm run playwright-setup` (browser binaries)
 *   - `pnpm run playwright-preflight-project` (Python env)
 *
 * Only run in an unsandboxed environment (requires localhost and a real browser).
 */

import test, { expect } from "@playwright/test";
import {
    createTemporaryProjectViaSyntheticAnndata,
    type SyntheticAnndataTemporaryProjectHandle,
} from "../utils/projectFixtures";

// Skipping the test so it doesn't run with all the other tests
test.describe.skip("Template", () => {
    test.setTimeout(180_000);

    test("example: project loads with a data count", async ({ page }) => {
        // createTemporaryProjectViaSyntheticAnndata creates a fresh project and
        // navigates `page` to it. Assign inside try; clean up in finally.
        let handle: SyntheticAnndataTemporaryProjectHandle | undefined;
        try {
            handle = await createTemporaryProjectViaSyntheticAnndata(page, {
                nameSegment: `template--${Date.now()}`,
                synthetic: {
                    profile: "minimal",
                    nCells: 200,
                    nGenes: 12,
                    force: true,
                },
            });

            // Your test logic here.
            // `page` is already on the project URL.
            const count = await page.locator(".ciview-cell-count").textContent();
            expect(count).toBeTruthy();
        } finally {
            if (handle !== undefined) {
                await handle.cleanup();
            }
        }
    });
});
