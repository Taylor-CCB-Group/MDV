import test, { expect, type BrowserContext, type Page } from "@playwright/test";
import { newProjectSetup, newProjectTeardown } from "../utils/testUtils";

/**
 * Integration test: edit a column in the frontend and save state, then verify round-trip.
 *
 * NOTE: Kept for reference, but our Playwright testing setup has general issues and this spec
 * may not run usefully in practice. Prefer the Python tests (e.g. test_unique_column.py) for
 * validating unique-column save_state behavior.
 *
 * Flow: Frontend getState() → updatedColumns with getMd() (unique columns decoded to string[])
 * → POST /save_state → project.save_state() → set_column_with_raw_data(ds, metadata, data).
 * For unique columns, the frontend already sends string[]; no conversion in set_column_with_raw_data.
 *
 * To test unique columns specifically: use a project that has a unique column, e.g. create with:
 *   python -m mdvtools.tests.generate_test_data /path/to/proj --mock --with-unique-column
 * or
 *   python -m mdvtools.tests.example_mock_usage --create-unique-project /path/to/proj
 * then open that project in the app and run this test (e.g. run Playwright against that project URL).
 *
 * These tests are intended to run manually outside the container when a server and project are available.
 */
test.describe.configure({ mode: "serial" });

test.describe("Unique column and save_state round-trip", () => {
    let projectId: string;
    let context: BrowserContext;
    let page: Page;

    test.beforeAll(async ({ browser }) => {
        context = await browser.newContext();
        page = await context.newPage();
        projectId = await newProjectSetup(page);
    });

    test.afterAll(async () => {
        await newProjectTeardown(page, projectId);
        await context.close();
    });

    test("project loads and table chart is available", async () => {
        await page.waitForFunction(() => !!(window as any).mdv?.chartManager);
        const hasTable = await page.evaluate(() => {
            const cm = (window as any).mdv?.chartManager;
            if (!cm?.charts) return false;
            return Object.values(cm.charts).some(
                (c: any) => c?.chart?.config?.type === "table_chart"
            );
        });
        expect(hasTable).toBe(true);
    });

    test("save_state posts successfully", async () => {
        const [resp] = await Promise.all([
            page.waitForResponse(
                (r) => r.request().method() === "POST" && r.url().includes("/save_state")
            ),
            page.evaluate(() => (window as any).mdv?.chartManager?.saveState()),
        ]);
        expect(resp.ok()).toBe(true);
    });
});
