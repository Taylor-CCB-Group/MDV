/**
 * Targets the Vite dev frontend (React color legend). Run with:
 *   TEST_BASE_URL=http://127.0.0.1:5170 pnpm run playwright-test-project tests_playwright/project/charts/color_legend_categorical.spec.ts --project=chromium --reporter=list
 */
import test, { expect } from "@playwright/test";
import {
    addChartViaUi,
    dragChartColorLegend,
    getChartColorLegendHost,
    getChartColorLegendState,
    isWorktreeFrontendForPlaywright,
    saveCurrentView,
    setChartColorBy,
    waitForChartByTitle,
} from "../../utils/helpers";
import {
    createTemporaryProjectViaSyntheticAnndata,
    waitForProjectReady,
    type SyntheticAnndataTemporaryProjectHandle,
} from "../../utils/projectFixtures";

test.describe("Color legend categorical", () => {
    test.setTimeout(180_000);

    test.beforeEach(() => {
        test.skip(
            !isWorktreeFrontendForPlaywright(),
            "React color legend requires the Vite dev server; set TEST_BASE_URL=http://127.0.0.1:5170",
        );
    });

    test("categorical legend visibility, labels, drag, and persistence after save reload", async ({
        page,
    }) => {
        const scatterTitle = `Scatter Color Legend ${Date.now()}`;

        let handle: SyntheticAnndataTemporaryProjectHandle | undefined;
        try {
            handle = await createTemporaryProjectViaSyntheticAnndata(page, {
                nameSegment: `color-legend-categorical--${Date.now()}`,
                synthetic: {
                    profile: "minimal",
                    nCells: 200,
                    nGenes: 12,
                    force: true,
                },
            });

            await addChartViaUi(page, "2D Scatter Plot", scatterTitle);
            await waitForChartByTitle(page, scatterTitle);
            await setChartColorBy(page, scatterTitle, "cell_type");

            const legend = getChartColorLegendHost(page, scatterTitle);
            await expect(legend).toBeVisible();

            const initial = await getChartColorLegendState(page, scatterTitle);
            expect(initial.columnTitle).toBe("cell_type");
            expect(initial.items.length).toBeGreaterThanOrEqual(2);

            const labels = initial.items.map((item) => item.label).sort();
            expect(labels).toEqual(expect.arrayContaining(["B-cell", "T-cell"]));
            for (const item of initial.items) {
                expect(item.color).toMatch(/^(#[0-9a-f]{3,8}|rgb\(|rgba\()/i);
            }

            const colorsByLabel = Object.fromEntries(initial.items.map((item) => [item.label, item.color]));

            const beforeDrag = { left: initial.left, top: initial.top };
            await dragChartColorLegend(page, scatterTitle, 100, 60);
            await expect
                .poll(async () => {
                    const moved = await getChartColorLegendState(page, scatterTitle);
                    return Math.abs(moved.left - beforeDrag.left) + Math.abs(moved.top - beforeDrag.top);
                })
                .toBeGreaterThan(20);

            await saveCurrentView(page);
            const saved = await getChartColorLegendState(page, scatterTitle);
            await page.goto(handle.projectUrl);
            await waitForProjectReady(page);

            await expect(getChartColorLegendHost(page, scatterTitle)).toBeVisible({ timeout: 60_000 });
            const reloaded = await getChartColorLegendState(page, scatterTitle);

            expect(reloaded.columnTitle).toBe("cell_type");
            expect(reloaded.items.map((item) => item.label).sort()).toEqual(
                saved.items.map((item) => item.label).sort(),
            );
            for (const item of reloaded.items) {
                expect(item.color).toBe(colorsByLabel[item.label]);
            }
            expect(Math.abs(reloaded.left - saved.left)).toBeLessThanOrEqual(3);
            expect(Math.abs(reloaded.top - saved.top)).toBeLessThanOrEqual(3);
        } finally {
            if (handle !== undefined) {
                await handle.cleanup();
            }
        }
    });
});
