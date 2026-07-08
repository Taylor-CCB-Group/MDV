/**
 * Targets the Vite dev frontend (React color legend). Run with:
 *   TEST_BASE_URL=http://127.0.0.1:5170 pnpm run playwright-test-project tests_playwright/project/charts/color_legend_continuous.spec.ts --project=chromium --reporter=list
 */
import test, { expect, type Page } from "@playwright/test";
import {
    addChartViaUi,
    dragChartColorLegend,
    getChartColorLegendHost,
    getChartColorLegendState,
    getFirstNumericColorByColumn,
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

async function getContinuousColorLegendFilterRuntimeState(page: Page, chartTitle: string) {
    return await page.evaluate((title) => {
        const chartManager = (window as any).mdv?.chartManager;
        const chartEntries = Object.values(chartManager?.charts ?? {});
        const chart = chartEntries
            .map((entry: any) => entry?.chart)
            .find((entry: any) => entry?.config?.title === title);
        if (!chart) {
            throw new Error(`Chart not found: ${title}`);
        }

        const panels = [...document.querySelectorAll(".ciview-chart-panel")];
        const panel = panels.find((entry) => {
            const titleNode = entry.querySelector(".ciview-chart-title");
            return titleNode?.textContent?.trim() === title;
        });

        return {
            chartFilter: chart.colorLegendFilter ?? null,
            dataStoreSize: chart.dataStore.size,
            filterSize: chart.dataStore.filterSize,
            resetDisplay: chart.resetButton?.style?.display ?? "",
            hasSelectedRange: (() => {
                const selection = panel?.querySelector(
                    ".legend-container .brush .selection",
                );
                if (!(selection instanceof SVGGraphicsElement)) {
                    return false;
                }
                const box = selection.getBBox();
                return getComputedStyle(selection).display !== "none" && box.width > 0;
            })(),
        };
    }, chartTitle);
}

async function dragContinuousColorLegendRangeOutsideRight(page: Page, chartTitle: string) {
    const slider = getChartColorLegendHost(page, chartTitle).locator(".brush .overlay");
    await expect(slider).toBeVisible();
    const box = await slider.boundingBox();
    if (!box) {
        throw new Error(`Continuous color legend slider not found for chart "${chartTitle}"`);
    }

    const startX = box.x + box.width * 0.35;
    const y = box.y + box.height / 2;
    await page.mouse.move(startX, y);
    await page.mouse.down();
    await page.mouse.move(box.x + box.width * 0.7, y);
    await page.mouse.move(box.x + box.width + 30, y);
    await page.mouse.up();
}

async function clearContinuousColorLegendRange(page: Page, chartTitle: string) {
    const slider = getChartColorLegendHost(page, chartTitle).locator(".brush .overlay");
    await expect(slider).toBeVisible();
    const box = await slider.boundingBox();
    if (!box) {
        throw new Error(`Continuous color legend slider not found for chart "${chartTitle}"`);
    }

    await page.mouse.click(box.x + box.width * 0.05, box.y + box.height / 2);
}

test.describe("Color legend continuous", () => {
    test.setTimeout(180_000);

    test.beforeEach(() => {
        test.skip(
            !isWorktreeFrontendForPlaywright(),
            "React color legend requires the Vite dev server; set TEST_BASE_URL=http://127.0.0.1:5170",
        );
    });

    test("continuous legend visibility, drag, and persistence after save reload", async ({
        page,
    }) => {
        const scatterTitle = `Scatter Continuous Color Legend ${Date.now()}`;

        let handle: SyntheticAnndataTemporaryProjectHandle | undefined;
        try {
            handle = await createTemporaryProjectViaSyntheticAnndata(page, {
                nameSegment: `color-legend-continuous--${Date.now()}`,
                synthetic: {
                    profile: "minimal",
                    nCells: 200,
                    nGenes: 12,
                    force: true,
                },
            });

            await addChartViaUi(page, "2D Scatter Plot", scatterTitle);
            await waitForChartByTitle(page, scatterTitle);

            const numericColorBy = await getFirstNumericColorByColumn(page);
            await setChartColorBy(page, scatterTitle, numericColorBy);

            const legend = getChartColorLegendHost(page, scatterTitle);
            await expect(legend).toBeVisible();
            await expect(legend.locator(".legend-body")).toHaveCount(0);

            const initial = await getChartColorLegendState(page, scatterTitle);
            expect(initial.width).toBeGreaterThan(40);
            expect(initial.height).toBeGreaterThan(20);

            const unfiltered = await getContinuousColorLegendFilterRuntimeState(
                page,
                scatterTitle,
            );
            expect(unfiltered.filterSize).toBe(unfiltered.dataStoreSize);
            expect(unfiltered.hasSelectedRange).toBe(false);

            await dragContinuousColorLegendRangeOutsideRight(page, scatterTitle);
            await expect
                .poll(async () => {
                    const state = await getContinuousColorLegendFilterRuntimeState(
                        page,
                        scatterTitle,
                    );
                    return (
                        state.chartFilter?.kind === "continuous" &&
                        state.filterSize < state.dataStoreSize &&
                        state.resetDisplay === "inline" &&
                        state.hasSelectedRange
                    );
                })
                .toBe(true);

            const filtered = await getContinuousColorLegendFilterRuntimeState(
                page,
                scatterTitle,
            );
            expect(filtered.filterSize).toBeLessThan(filtered.dataStoreSize);

            await clearContinuousColorLegendRange(page, scatterTitle);
            await expect
                .poll(async () => {
                    const state = await getContinuousColorLegendFilterRuntimeState(
                        page,
                        scatterTitle,
                    );
                    return (
                        state.chartFilter === null &&
                        state.filterSize === unfiltered.dataStoreSize &&
                        state.resetDisplay === "none" &&
                        !state.hasSelectedRange
                    );
                })
                .toBe(true);

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
            expect(reloaded.width).toBeGreaterThan(40);
            expect(reloaded.height).toBeGreaterThan(20);
            expect(Math.abs(reloaded.left - saved.left)).toBeLessThanOrEqual(3);
            expect(Math.abs(reloaded.top - saved.top)).toBeLessThanOrEqual(3);
        } finally {
            if (handle !== undefined) {
                await handle.cleanup();
            }
        }
    });
});
