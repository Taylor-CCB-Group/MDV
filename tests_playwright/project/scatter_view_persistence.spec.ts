import test, { expect } from "@playwright/test";
import { createTemporaryProject, waitForProjectReady } from "../utils/projectFixtures";
import {
    addChartViaUi,
    GENERATED_MOCK_PROJECT,
    getChartSummaries,
    saveCurrentView,
    waitForChartByTitle,
    waitForViewUnsavedState,
} from "./helpers";

test.describe("Scatter View Persistence", () => {
    test.setTimeout(180_000);
    test("creates a 2D Scatter Plot, saves the view, and keeps it after reload", async ({ page }) => {
        const projectHandle = await createTemporaryProject(page, {
            allowCsvFallback: false,
            mockConfig: GENERATED_MOCK_PROJECT,
            projectName: "Scatter View Persistence Fixture",
        });

        try {
            const chartsBefore = await getChartSummaries(page);
            const chartTitle = `Playwright Scatter ${Date.now()}`;

            await addChartViaUi(page, "2D Scatter Plot", chartTitle);
            await waitForViewUnsavedState(page, true);
            await waitForChartByTitle(page, chartTitle);

            const chartsAfterCreate = await getChartSummaries(page);
            const createdChart = chartsAfterCreate.find(
                (chart) => chart.title === chartTitle && !chartsBefore.some((before) => before.id === chart.id),
            );

            expect(createdChart).toBeTruthy();
            expect(createdChart?.type.toLowerCase()).toContain("scatter");

            await saveCurrentView(page);

            await page.reload();
            await waitForProjectReady(page);

            const chartsAfterReload = await getChartSummaries(page);
            expect(
                chartsAfterReload.some(
                    (chart) =>
                        chart.id === createdChart?.id &&
                        chart.title === chartTitle &&
                        chart.type === createdChart.type,
                ),
            ).toBe(true);
        } finally {
            await projectHandle.cleanup();
        }
    });
});
