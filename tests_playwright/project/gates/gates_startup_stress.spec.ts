import test, { expect } from "@playwright/test";
import {
    createTemporaryProjectViaSyntheticAnndata,
    type SyntheticAnndataTemporaryProjectHandle,
} from "../../utils/projectFixtures";

// Pass the PLAYWRIGHT_STRESS=1 to run this test
test.describe("gates startup stress", () => {
    test.skip(
        process.env.PLAYWRIGHT_STRESS !== "1",
        "10M-row startup stress test is opt-in. Set PLAYWRIGHT_STRESS=1 to run.",
    );

    test.setTimeout(600_000);

    test("10M-cell project stays alive during startup", async ({ page }) => {
        let handle: SyntheticAnndataTemporaryProjectHandle | undefined;
        const pageErrors: string[] = [];
        const consoleErrors: string[] = [];

        page.on("pageerror", (error) => {
            pageErrors.push(error.message);
        });
        page.on("console", (message) => {
            if (message.type() === "error") {
                consoleErrors.push(message.text());
            }
        });

        try {
            handle = await createTemporaryProjectViaSyntheticAnndata(page, {
                nameSegment: `synth-anndata--gates-startup-stress--10m--g1--${Date.now()}`,
                synthetic: {
                    profile: "memory-efficient",
                    nCells: 10_000_000,
                    nGenes: 1,
                    sparseDensity: 0.01,
                    force: true,
                },
            });

            await expect(page.locator(".ciview-contentDiv").first()).toBeVisible();

            await page.waitForTimeout(60_000);

            await expect(page.locator(".ciview-contentDiv").first()).toBeVisible();
            const gatesColumnState = await page.evaluate(() => {
                const mdv = Reflect.get(window, "mdv");
                const chartManager = mdv?.chartManager;
                const dataStore =
                    chartManager?.getDataSource?.("cells") ??
                    chartManager?.dsIndex?.cells?.dataStore;
                const gatesColumn = dataStore?.columnIndex?.__gates__;

                return gatesColumn
                    ? {
                        hasColumn: true,
                        hasData: Boolean(gatesColumn.data),
                        datatype: gatesColumn.datatype,
                    }
                    : { hasColumn: false };
            });

            expect(page.url()).toContain(`/project/${handle.projectId}`);
            expect(await page.title()).not.toContain("crashed");
            expect(gatesColumnState).toEqual({ hasColumn: false });
            expect(pageErrors).toEqual([]);
            expect(consoleErrors.filter((text) => text.includes("__gates__"))).toEqual([]);
        } finally {
            if (handle !== undefined) {
                await handle.cleanup();
            }
        }
    });
});
