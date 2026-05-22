import test, { expect } from "@playwright/test";
import { createTemporaryProjectViaSyntheticAnndata } from "../../utils/projectFixtures";
import {
    addChartViaUi,
    closeTopDialog,
    expectChartCreationErrorVisible,
    expectChartPanelHasNoError,
    getChartSummaries,
    openChartDebugDialog,
    openChartSettingsDialog,
    waitForChartByTitle,
    waitForViewUnsavedState,
} from "../../utils/helpers";

test.describe("Chart Creation Single", () => {
    test.setTimeout(180_000);

    test("creates one 2D scatter chart smoke path using synthetic AnnData project", async ({
        page,
    }) => {
        const projectHandle = await createTemporaryProjectViaSyntheticAnndata(page, {
            synthetic: {
                profile: "minimal",
                nCells: 200,
                nGenes: 12,
                force: true,
            },
        });

        try {
            expect(projectHandle.sourceUsed).toBe("synthetic-anndata-filesystem-rescan");

            const title = `2D Scatter Plot ${Date.now()}`;
            const chartsBefore = await getChartSummaries(page);

            await addChartViaUi(page, "2D Scatter Plot", title);
            await waitForViewUnsavedState(page, true);
            await waitForChartByTitle(page, title);
            await expectChartPanelHasNoError(page, title);

            const chartsAfter = await getChartSummaries(page);
            const createdChart = chartsAfter.find(
                (chart) =>
                    chart.title === title &&
                    ["wgl_scatter_plot", "wgl_scatter_plot_dev"].includes(chart.type) &&
                    !chartsBefore.some((before) => before.id === chart.id),
            );

            expect(createdChart).toBeTruthy();
            expect(createdChart?.type).toBeTruthy();

            await openChartSettingsDialog(page, title);
            await expect(page.getByLabel("Search Settings by Folder or Name")).toBeVisible();
            await closeTopDialog(page);

            await openChartDebugDialog(page, title);
            await expect(page.getByPlaceholder("Filter...")).toBeVisible();
            await closeTopDialog(page);
        } finally {
            await projectHandle.cleanup();
        }
    });

    test("shows a chart error overlay when chart creation fails", async ({ page }) => {
        const projectHandle = await createTemporaryProjectViaSyntheticAnndata(page, {
            synthetic: {
                profile: "minimal",
                nCells: 200,
                nGenes: 12,
                force: true,
            },
        });

        try {
            const title = `2D Scatter Plot Forced Error ${Date.now()}`;

            await page.evaluate(() => {
                const chartManager = (window as any).mdv?.chartManager;
                const originalAddChart = chartManager?._addChart;
                if (!chartManager || typeof originalAddChart !== "function") {
                    throw new Error("ChartManager _addChart hook is unavailable.");
                }

                chartManager._addChart = async function patchedAddChart(
                    ...args: unknown[]
                ) {
                    chartManager._addChart = originalAddChart;
                    throw new Error("Playwright forced chart creation failure");
                };
            });

            await addChartViaUi(page, "2D Scatter Plot", title);
            await expectChartCreationErrorVisible(page, "Playwright forced chart creation failure");
        } finally {
            await projectHandle.cleanup();
        }
    });
});
