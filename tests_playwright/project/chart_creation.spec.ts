import test, { expect } from "@playwright/test";
import {
    createSharedSyntheticAnndataSuite,
    type SharedSyntheticAnndataSuiteHandle,
} from "../utils/projectFixtures";
import {
    addChartViaUi,
    getChartSummaries,
    waitForChartByTitle,
    waitForViewUnsavedState,
} from "../utils/helpers";

type ChartCreationCase = {
    chartName: string;
    /** Acceptable `chart.type` values after creation (synthetic AnnData may vary for scatter). */
    expectedTypes: readonly string[];
};

const CHART_CASES: ChartCreationCase[] = [
    {
        chartName: "2D Scatter Plot",
        expectedTypes: ["wgl_scatter_plot", "wgl_scatter_plot_dev"],
    },
    {
        chartName: "Row Chart",
        expectedTypes: ["row_chart"],
    },
    {
        chartName: "Table",
        expectedTypes: ["table_chart_react"],
    },
];

test.describe("Chart Creation", () => {
    test.describe.configure({ mode: "serial" });
    test.setTimeout(180_000);

    let sharedProjectHandle: SharedSyntheticAnndataSuiteHandle | undefined;

    test.beforeAll(async ({ browser }) => {
        sharedProjectHandle = await createSharedSyntheticAnndataSuite(browser, {
            nameSegment: `chart-creation--shared--synth-anndata--${Date.now()}`,
            synthetic: {
                profile: "minimal",
                nCells: 200,
                nGenes: 12,
                force: true,
            },
        });
        expect(sharedProjectHandle.sourceUsed).toBe("synthetic-anndata-filesystem-rescan");
    });

    test.afterAll(async () => {
        await sharedProjectHandle?.cleanup();
    });

    for (const chartCase of CHART_CASES) {
        test(`creates ${chartCase.chartName} in a shared synthetic AnnData project`, async ({ page }) => {
            expect(sharedProjectHandle).toBeTruthy();
            if (!sharedProjectHandle) {
                throw new Error("Shared synthetic AnnData project was not initialized.");
            }

            await sharedProjectHandle.openProjectPage(page);

            const title = `${chartCase.chartName} ${Date.now()}`;
            const chartsBefore = await getChartSummaries(page);

            await addChartViaUi(page, chartCase.chartName, title);
            await waitForViewUnsavedState(page, true);
            await waitForChartByTitle(page, title);

            const chartsAfter = await getChartSummaries(page);
            const createdChart = chartsAfter.find(
                (chart) =>
                    chart.title === title &&
                    chartCase.expectedTypes.includes(chart.type) &&
                    !chartsBefore.some((before) => before.id === chart.id),
            );

            expect(createdChart).toBeTruthy();
        });
    }
});
