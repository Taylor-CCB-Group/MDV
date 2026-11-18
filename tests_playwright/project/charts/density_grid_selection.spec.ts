/**
 * The deck scatter plot's density grid: one cell per density field, all sharing one view.
 * - each field gets a labelled cell, drawn without errors;
 * - a selection drawn in any cell lands where it would in the first, since every cell shows
 *   the same view (a regression guard for unprojecting through the wrong cell's viewport);
 *   the cells compared are in different rows and columns;
 * - the toolbar switches between the grid and the single plot.
 *
 * Runs against whichever frontend TEST_BASE_URL serves; to test a branch's UI, use Vite from that checkout, e.g.
 *   PORT=5171 pnpm run dev
 *   PLAYWRIGHT_HTML_OPEN=never TEST_BASE_URL=http://127.0.0.1:5171 pnpm run playwright-test-project tests_playwright/project/charts/density_grid_selection.spec.ts --project=chromium --reporter=list
 */

import test, { expect, type Locator, type Page } from "@playwright/test";
import { expectChartPanelHasNoError, waitForChartByTitle } from "../../utils/helpers";
import { createTemporaryProjectViaSyntheticAnndata } from "../../utils/projectFixtures";

type Point = [number, number];
type Box = { x: number; y: number; width: number; height: number };

const CHART_ATTR = "data-density-grid-test";
const FIELD_COUNT = 4;

/** A 2D scatter plot set to show `FIELD_COUNT` numeric columns as a density grid. */
async function addDensityGridChart(page: Page, title: string) {
    // addChart throws in positionChart until the view has loaded
    await page.waitForFunction(
        () => {
            const cm = (window as any).mdv?.chartManager;
            return Boolean(cm?.viewData?.dataSources) && !cm._inInit;
        },
        undefined,
        { timeout: 60_000 },
    );
    const id = `density-grid-${Date.now()}`;
    const fields: string[] = await page.evaluate(
        ({ id, title, fieldCount }) => {
            const cm = (window as any).mdv.chartManager;
            const dsName = Object.keys(cm.dsIndex)[0];
            const columns: { field: string; datatype: string }[] = cm.dsIndex[dsName].dataStore.columns;
            const numeric = columns
                .filter((c) => ["double", "integer", "int32"].includes(c.datatype))
                .map((c) => c.field);
            // density fields must be `double` (see `isLoadedDoubleColumn`)
            const doubles = columns.filter((c) => c.datatype === "double").map((c) => c.field);
            const umap = ["X_umap_1", "X_umap_2"];
            const param = umap.every((f) => numeric.includes(f)) ? umap : numeric.slice(0, 2);
            const densityFields = doubles.filter((f) => !param.includes(f)).slice(0, fieldCount);
            if (param.length < 2 || densityFields.length < fieldCount) {
                throw new Error(`need 2 numeric and ${fieldCount} more double columns in ${dsName}, got ${numeric}`);
            }
            cm.addChart(dsName, {
                type: "DeckContourScatter",
                id,
                title,
                param,
                densityFields,
                density_mode: "grid",
                size: [700, 700],
                position: [10, 10],
            });
            return densityFields;
        },
        { id, title, fieldCount: FIELD_COUNT },
    );
    await waitForChartByTitle(page, title);
    await page.waitForFunction(
        ({ id, attr }) => {
            const chart = (window as any).mdv.chartManager.charts[id]?.chart;
            if (!chart?.contentDiv) return false;
            chart.contentDiv.setAttribute(attr, id);
            return true;
        },
        { id, attr: CHART_ATTR },
        { timeout: 60_000 },
    );
    return { id, fields };
}

function chartScope(page: Page, id: string): Locator {
    return page.locator(`[${CHART_ATTR}="${id}"]`);
}

async function cellBox(page: Page, id: string, index: number): Promise<Box> {
    const cell = chartScope(page, id).locator(".mdv-chart-array__cell").nth(index);
    await cell.scrollIntoViewIfNeeded();
    const box = await cell.boundingBox();
    if (!box) throw new Error(`no box for cell ${index}`);
    expect(box.width).toBeGreaterThan(50);
    expect(box.height).toBeGreaterThan(50);
    return box;
}

function at(box: Box, fx: number, fy: number): Point {
    return [box.x + box.width * fx, box.y + box.height * fy];
}

async function readRings(page: Page, id: string): Promise<Point[][]> {
    return await page.evaluate((id) => {
        const config = (window as any).mdv.chartManager.charts[id].chart.config;
        const features: { geometry?: { type?: string; coordinates?: number[][][] } }[] =
            config.selectionFeatureCollection?.features ?? [];
        return features.map((f) =>
            f.geometry?.type === "Polygon" && f.geometry.coordinates
                ? f.geometry.coordinates[0].map((p): [number, number] => [p[0], p[1]])
                : [],
        );
    }, id);
}

async function clearSelection(page: Page, id: string) {
    await page.evaluate((id) => (window as any).mdv.chartManager.charts[id].chart.removeFilter(), id);
    await expect.poll(async () => (await readRings(page, id)).length).toBe(0);
}

/** Drags a rectangle over the middle of a cell and returns the ring it made, in data units. */
async function dragRectangleInCell(page: Page, id: string, cell: number): Promise<Point[]> {
    await clearSelection(page, id);
    const box = await cellBox(page, id, cell);
    await page.mouse.move(...at(box, 0.35, 0.35), { steps: 3 });
    await page.mouse.down();
    await page.mouse.move(...at(box, 0.65, 0.65), { steps: 10 });
    await page.mouse.up();
    await expect.poll(async () => (await readRings(page, id)).length, { timeout: 20_000 }).toBe(1);
    return (await readRings(page, id))[0];
}

function bounds(ring: Point[]) {
    const xs = ring.map((p) => p[0]);
    const ys = ring.map((p) => p[1]);
    return { minX: Math.min(...xs), maxX: Math.max(...xs), minY: Math.min(...ys), maxY: Math.max(...ys) };
}

test.describe("Density grid", () => {
    test.setTimeout(180_000);

    test("draws a cell per field, selects through any cell, and toggles back to one plot", async ({ page }) => {
        const projectHandle = await createTemporaryProjectViaSyntheticAnndata(page, {
            // "realistic" has several numeric obs columns and the X_umap / X_pca reductions
            synthetic: { profile: "realistic", nCells: 500, nGenes: 50, force: true },
        });

        try {
            const title = `Density Grid ${Date.now()}`;
            const { id, fields } = await addDensityGridChart(page, title);
            await expectChartPanelHasNoError(page, title);

            const scope = chartScope(page, id);
            const cells = scope.locator(".mdv-chart-array__cell");
            await expect(cells).toHaveCount(fields.length, { timeout: 30_000 });
            for (const field of fields) {
                await expect(scope.locator(".mdv-chart-array__cell", { hasText: field })).toHaveCount(1);
            }
            // a canvas of its own in every cell, drawn at the cell's size
            const cellCanvases = scope.locator(".mdv-chart-array__cell > canvas");
            await expect(cellCanvases).toHaveCount(fields.length);
            await expect(cellCanvases.first()).toBeVisible({ timeout: 30_000 });
            await expect
                .poll(async () => await cellCanvases.first().evaluate((canvas: HTMLCanvasElement) => canvas.width))
                .toBeGreaterThan(50);

            // the toolbar starts collapsed to the active tool
            await scope.locator('button[aria-label="Show tools"]').click();
            await scope.locator('button[aria-label="Rectangle"]').click();
            // the last cell is in another row and column (a 2 x 2 grid at this size)
            const last = fields.length - 1;
            expect((await cellBox(page, id, last)).y).toBeGreaterThan((await cellBox(page, id, 0)).y);
            const first = bounds(await dragRectangleInCell(page, id, 0));
            const second = bounds(await dragRectangleInCell(page, id, last));
            expect(first.maxX).toBeGreaterThan(first.minX);
            expect(first.maxY).toBeGreaterThan(first.minY);
            // the same place in another cell is the same place in the data; allow a pixel or two
            const tolerance = (first.maxX - first.minX) * 0.05;
            for (const key of ["minX", "maxX", "minY", "maxY"] as const) {
                expect(Math.abs(second[key] - first[key])).toBeLessThan(tolerance);
            }

            await scope.locator('button[aria-label="Density overlay"]').click();
            await expect(cells).toHaveCount(0);
            await expect(scope.getByLabel("Scatter plot")).toBeVisible();
            await scope.locator('button[aria-label="Density grid"]').click();
            await expect(cells).toHaveCount(fields.length);
        } finally {
            await projectHandle.cleanup();
        }
    });
});
