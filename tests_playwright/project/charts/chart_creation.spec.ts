import test, { expect } from "@playwright/test";
import {
    createSharedSyntheticAnndataSuite,
    type SharedSyntheticAnndataSuiteHandle,
} from "../../utils/projectFixtures";
import {
    addChartViaUi,
    getAvailableAddableChartTypes,
    getChartSummaries,
    waitForChartByTitle,
    waitForViewUnsavedState,
} from "../../utils/helpers";

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

    test("creates every addable chart in a shared synthetic AnnData project", async ({ page }) => {
        expect(sharedProjectHandle).toBeTruthy();
        if (!sharedProjectHandle) {
            throw new Error("Shared synthetic AnnData project was not initialized.");
        }

        await sharedProjectHandle.openProjectPage(page);

        const availableChartTypes = await getAvailableAddableChartTypes(page);
        expect(availableChartTypes.length).toBeGreaterThan(0);

        for (const [index, chartCase] of availableChartTypes.entries()) {
            await test.step(`creates ${chartCase.chartName}`, async () => {
                const title = `${chartCase.chartName} ${index + 1} ${Date.now()}`;
                const chartsBefore = await getChartSummaries(page);

                await addChartViaUi(page, chartCase.chartName, title);
                await waitForViewUnsavedState(page, true);
                await waitForChartByTitle(page, title);

                const chartsAfter = await getChartSummaries(page);
                const createdChart = chartsAfter.find(
                    (chart) =>
                        chart.title === title &&
                        !chartsBefore.some((before) => before.id === chart.id),
                );

                expect(createdChart).toBeTruthy();
                expect(createdChart?.type).toBeTruthy();
            });
        }
    });
});
