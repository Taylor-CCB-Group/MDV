import test, { expect } from "@playwright/test";
import { createTemporaryProjectViaSyntheticAnndata } from "../utils/projectFixtures";
import {
    addChartViaUi,
    getChartSummaries,
    waitForChartByTitle,
    waitForViewUnsavedState,
} from "../utils/helpers";

test.describe("Chart Creation Single", () => {
    test.setTimeout(180_000);

    test("creates one 2D scatter chart using synthetic AnnData project (filesystem rescan)", async ({
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

            const chartsAfter = await getChartSummaries(page);
            const createdChart = chartsAfter.find(
                (chart) =>
                    chart.title === title &&
                    ["wgl_scatter_plot", "wgl_scatter_plot_dev"].includes(chart.type) &&
                    !chartsBefore.some((before) => before.id === chart.id),
            );

            expect(createdChart).toBeTruthy();
        } finally {
            await projectHandle.cleanup();
        }
    });
});
