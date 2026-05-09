/**
 * One synthetic AnnData project (filesystem + `/rescan_projects`) shared across tests so
 * Python generation and rescan run once; tests run serially and add charts in order.
 *
 * The first test uses the Playwright `page` fixture to create the project (same as
 * `chart_creation_single.spec.ts`), so `baseURL` from `playwright.config.ts` applies to
 * `page.request` and `goto("/project/...")`. A `beforeAll` + `browser.newPage()` pattern
 * does *not* get that `baseURL` unless you use `browser.newContext({ baseURL })`.
 */
import test, { expect } from "@playwright/test";
import {
    createTemporaryProjectViaSyntheticAnndata,
    type SyntheticAnndataTemporaryProjectHandle,
    waitForProjectReady,
} from "../utils/projectFixtures";
import {
    addChartViaUi,
    getChartSummaries,
    waitForChartByTitle,
    waitForViewUnsavedState,
} from "./helpers";

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

    test.setTimeout(300_000);

    let projectHandle: SyntheticAnndataTemporaryProjectHandle | undefined;

    test.afterAll(async () => {
        if (projectHandle) {
            await projectHandle.cleanup();
        }
    });

    for (const chartCase of CHART_CASES) {
        test(`creates ${chartCase.chartName} in shared synthetic AnnData project`, async ({ page }) => {
            if (!projectHandle) {
                projectHandle = await createTemporaryProjectViaSyntheticAnndata(page, {
                    nameSegment: `chart-creation--synth-anndata--shared--${Date.now()}`,
                    synthetic: {
                        profile: "minimal",
                        nCells: 200,
                        nGenes: 12,
                        force: true,
                    },
                });
            } else {
                await page.goto(projectHandle.projectUrl);
                await waitForProjectReady(page);
            }

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
