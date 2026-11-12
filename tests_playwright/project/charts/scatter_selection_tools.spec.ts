/**
 * End-to-end guard for the deck scatter selection tools (issues #577, #578):
 * - a double-click finishes a Polygon with the double-clicked point as its last vertex, and doesn't zoom;
 * - Freehand, Modify (vertex removal) and Pan (translate selection / controller pan + zoom) work;
 * - all of that still works after the chart is popped out into another window, and after it comes back.
 *
 * Runs against whichever frontend TEST_BASE_URL serves; to test a branch's UI, use Vite from that checkout, e.g.
 *   PORT=5171 pnpm run dev
 *   PLAYWRIGHT_HTML_OPEN=never TEST_BASE_URL=http://127.0.0.1:5171 pnpm run playwright-test-project tests_playwright/project/charts/scatter_selection_tools.spec.ts --project=chromium --reporter=list
 */

import test, { expect, type Locator, type Page } from "@playwright/test";
import { waitForChartByTitle } from "../../utils/helpers";
import {
    createTemporaryProjectViaSyntheticAnndata,
    type SyntheticAnndataTemporaryProjectHandle,
} from "../../utils/projectFixtures";

type Point = [number, number];
type Box = { x: number; y: number; width: number; height: number };
type SelectionViewState = { target: number[]; zoom: number };
type SelectionState = { rings: Point[][]; viewState: SelectionViewState | null };
type Tool = "Pan" | "Rectangle" | "Polygon" | "Freehand" | "Modify";

const CHART_ATTR = "data-selection-tools-test";
/** Gap between separate clicks, so they aren't read as a double-click and picking can keep up. */
const CLICK_GAP_MS = 650;

async function addScatterChart(page: Page, title: string): Promise<string> {
    // addChart throws in positionChart until the view has loaded (its view may have no charts to wait for)
    await page.waitForFunction(
        () => {
            const cm = (window as any).mdv?.chartManager;
            return Boolean(cm?.viewData?.dataSources) && !cm._inInit;
        },
        undefined,
        { timeout: 60_000 },
    );
    const id = `selection-tools-${Date.now()}`;
    await page.evaluate(
        ({ id, title }) => {
            const cm = (window as any).mdv.chartManager;
            const dsName = Object.keys(cm.dsIndex)[0];
            const columns: { field: string; datatype: string }[] = cm.dsIndex[dsName].dataStore.columns;
            const numeric = columns
                .filter((c) => ["double", "integer", "int32"].includes(c.datatype))
                .map((c) => c.field);
            const umap = ["X_umap_1", "X_umap_2"];
            const param = umap.every((f) => numeric.includes(f)) ? umap : numeric.slice(0, 2);
            if (param.length < 2) throw new Error(`need two numeric columns in ${dsName}, got ${numeric}`);
            cm.addChart(dsName, {
                // DECK_CONTOUR_SCATTER_TYPE: what Add Chart's "2D Scatter Plot" creates
                type: "DeckContourScatter",
                id,
                title,
                param,
                size: [700, 700],
                position: [10, 10],
            });
        },
        { id, title },
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
    return id;
}

function chartScope(ui: Page, id: string): Locator {
    return ui.locator(`[${CHART_ATTR}="${id}"]`);
}

async function selectTool(ui: Page, id: string, tool: Tool) {
    await chartScope(ui, id).locator(`button[aria-label="${tool}"]`).click();
}

/** The deck canvas is the largest canvas in the chart's content div. */
async function measureCanvas(ui: Page, id: string): Promise<Box> {
    const canvases = chartScope(ui, id).locator("canvas");
    await expect.poll(() => canvases.count()).toBeGreaterThan(0);
    const box = await canvases.evaluateAll((els) => {
        const rects = els.map((el) => el.getBoundingClientRect());
        const r = rects.reduce((a, b) => (b.width * b.height > a.width * a.height ? b : a));
        return { x: r.x, y: r.y, width: r.width, height: r.height };
    });
    expect(box.width).toBeGreaterThan(100);
    expect(box.height).toBeGreaterThan(100);
    return box;
}

/** Selection and view state live on the chart config in the main window, wherever the chart is shown. */
async function readState(page: Page, id: string): Promise<SelectionState> {
    return await page.evaluate((id) => {
        const config = (window as any).mdv.chartManager.charts[id].chart.config;
        const features: { geometry?: { type?: string; coordinates?: number[][][] } }[] =
            config.selectionFeatureCollection?.features ?? [];
        const rings = features.map((f) =>
            f.geometry?.type === "Polygon" && f.geometry.coordinates
                ? f.geometry.coordinates[0].map((p): [number, number] => [p[0], p[1]])
                : [],
        );
        const vs = config.viewState;
        return {
            rings,
            viewState: vs ? { target: [...vs.target], zoom: vs.zoom } : null,
        };
    }, id);
}

async function clearSelection(page: Page, id: string) {
    await page.evaluate((id) => (window as any).mdv.chartManager.charts[id].chart.removeFilter(), id);
    await expect.poll(async () => (await readState(page, id)).rings.length).toBe(0);
}

async function waitForViewState(page: Page, id: string): Promise<SelectionViewState> {
    await expect.poll(async () => (await readState(page, id)).viewState).not.toBeNull();
    const { viewState } = await readState(page, id);
    if (!viewState) throw new Error("no viewState");
    return viewState;
}

/** A point in the canvas given as fractions of its width / height. */
function at(box: Box, fx: number, fy: number): Point {
    return [box.x + box.width * fx, box.y + box.height * fy];
}

/** Five distinct, well separated vertices; the last is the double-clicked one. */
function pentagon(box: Box): Point[] {
    return Array.from({ length: 5 }, (_, i) => {
        const a = -Math.PI / 2 + (i * 2 * Math.PI) / 5;
        return at(box, 0.5 + 0.25 * Math.cos(a), 0.5 + 0.25 * Math.sin(a));
    });
}

async function pause(ui: Page, ms = CLICK_GAP_MS) {
    await ui.waitForTimeout(ms);
}

/** Software-rendered picking is slow on the first interaction; get that out of the way before asserting. */
async function warmUp(ui: Page, page: Page, id: string, box: Box) {
    await selectTool(ui, id, "Freehand");
    await ui.mouse.move(...at(box, 0.4, 0.4), { steps: 3 });
    await ui.mouse.down();
    for (const p of [at(box, 0.6, 0.4), at(box, 0.6, 0.6), at(box, 0.4, 0.6), at(box, 0.4, 0.4)]) {
        await ui.mouse.move(...p, { steps: 5 });
    }
    await ui.mouse.up();
    await expect.poll(async () => (await readState(page, id)).rings.length, { timeout: 20_000 }).toBe(1);
    await clearSelection(page, id);
}

/** (a) + (b): returns the clicked screen points and the resulting ring for calibration. */
async function drawPolygonWithDoubleClick(ui: Page, page: Page, id: string, box: Box) {
    await clearSelection(page, id);
    await selectTool(ui, id, "Polygon");
    const before = await waitForViewState(page, id);
    const points = pentagon(box);
    for (const [x, y] of points.slice(0, 4)) {
        await ui.mouse.move(x, y, { steps: 3 });
        await ui.mouse.click(x, y);
        await pause(ui);
    }
    const [lx, ly] = points[4];
    await ui.mouse.move(lx, ly, { steps: 3 });
    await ui.mouse.dblclick(lx, ly);
    await expect.poll(async () => (await readState(page, id)).rings.map((r) => r.length)).toEqual([6]);
    const { rings, viewState } = await readState(page, id);
    // (b) doubleClickZoom is off: the double-click finishes the polygon and leaves the view alone
    expect(viewState?.zoom).toBe(before.zoom);
    expect(viewState?.target).toEqual(before.target);
    return { points, ring: rings[0] };
}

/** (c) */
async function drawFreehandSquare(ui: Page, page: Page, id: string, box: Box) {
    await clearSelection(page, id);
    await selectTool(ui, id, "Freehand");
    const corners = [at(box, 0.3, 0.3), at(box, 0.7, 0.3), at(box, 0.7, 0.7), at(box, 0.3, 0.7)];
    await ui.mouse.move(...corners[0], { steps: 3 });
    await ui.mouse.down();
    for (const corner of [...corners.slice(1), corners[0]]) {
        await ui.mouse.move(...corner, { steps: 8 });
        await pause(ui, 150);
    }
    await ui.mouse.up();
    await expect.poll(async () => (await readState(page, id)).rings.length).toBe(1);
    const { rings } = await readState(page, id);
    expect(rings[0].length).toBeGreaterThan(4);
    return rings[0];
}

/** (f): with no selection, Pan drags the view and the wheel zooms. */
async function checkController(ui: Page, page: Page, id: string, box: Box) {
    await clearSelection(page, id);
    await selectTool(ui, id, "Pan");
    const before = await waitForViewState(page, id);
    const [sx, sy] = at(box, 0.4, 0.4);
    await ui.mouse.move(sx, sy, { steps: 3 });
    await ui.mouse.down();
    await ui.mouse.move(sx + 120, sy + 80, { steps: 10 });
    await ui.mouse.up();
    await expect.poll(async () => (await readState(page, id)).viewState?.target).not.toEqual(before.target);
    const panned = await waitForViewState(page, id);
    expect(panned.zoom).toBe(before.zoom);

    const [cx, cy] = at(box, 0.5, 0.5);
    await ui.mouse.move(cx, cy, { steps: 3 });
    await ui.mouse.wheel(0, -400);
    await expect.poll(async () => (await readState(page, id)).viewState?.zoom).not.toBe(before.zoom);
}

type DataToScreen = (p: Point) => Point;

/** (d): Modify removes a clicked vertex of the selection. Returns the edited ring. */
async function removeVertex(ui: Page, page: Page, id: string, box: Box, toScreen: DataToScreen, ring: Point[]) {
    await selectTool(ui, id, "Modify");
    await ui.mouse.move(...at(box, 0.5, 0.5), { steps: 5 });
    await pause(ui);
    // a vertex away from where the freehand started and finished
    const vertex = toScreen(ring[Math.floor(ring.length / 2)]);
    await ui.mouse.move(...vertex, { steps: 5 });
    await pause(ui);
    await ui.mouse.click(...vertex);
    await expect.poll(async () => (await readState(page, id)).rings[0]?.length).toBe(ring.length - 1);
    return (await readState(page, id)).rings[0];
}

/** (e): with a selection, Pan drags it - every vertex by the same delta - and leaves the view alone. */
async function translateSelection(ui: Page, page: Page, id: string, box: Box, ring: Point[]) {
    await selectTool(ui, id, "Pan");
    const view = await waitForViewState(page, id);
    const [sx, sy] = at(box, 0.5, 0.5);
    await ui.mouse.move(sx, sy, { steps: 5 });
    await pause(ui);
    await ui.mouse.down();
    await ui.mouse.move(sx + 60, sy + 40, { steps: 10 });
    await ui.mouse.up();
    await expect.poll(async () => (await readState(page, id)).rings[0]?.[0]).not.toEqual(ring[0]);
    const { rings, viewState } = await readState(page, id);
    expect(rings[0].length).toBe(ring.length);
    const [dx, dy] = [rings[0][0][0] - ring[0][0], rings[0][0][1] - ring[0][1]];
    expect(Math.hypot(dx, dy)).toBeGreaterThan(0);
    for (let k = 0; k < ring.length; k++) {
        expect(rings[0][k][0] - ring[k][0]).toBeCloseTo(dx, 6);
        expect(rings[0][k][1] - ring[k][1]).toBeCloseTo(dy, 6);
    }
    expect(viewState?.target).toEqual(view.target);
}

/**
 * (a) to (f) in whichever window the chart is in: polygon (which also calibrates data to screen),
 * freehand, modify, drag the selection, then pan and zoom the view.
 */
async function checkAllTools(ui: Page, page: Page, id: string, box: Box, where: string) {
    let toScreen: DataToScreen | undefined;
    await test.step(`${where}: (a, b) polygon finishes on double-click, keeping its point, without zooming`, async () => {
        const { points, ring } = await drawPolygonWithDoubleClick(ui, page, id, box);
        expect(ring[4][0]).not.toBe(ring[3][0]);
        toScreen = makeDataToScreen(points, ring);
    });
    let square: Point[] = [];
    await test.step(`${where}: (c) freehand drag draws one polygon`, async () => {
        square = await drawFreehandSquare(ui, page, id, box);
    });
    await test.step(`${where}: (d) modify removes the clicked vertex`, async () => {
        if (!toScreen) throw new Error("no calibration from (a)");
        square = await removeVertex(ui, page, id, box, toScreen, square);
    });
    await test.step(`${where}: (e) pan with a selection translates every vertex by the same delta`, async () => {
        await translateSelection(ui, page, id, box, square);
    });
    await test.step(`${where}: (f) pan and wheel move the view`, async () => {
        await checkController(ui, page, id, box);
    });
}

/**
 * Map data coordinates to screen using the polygon drawn in (a): clicks and vertices correspond in
 * order, and the view is orthographic so each axis is linear.
 */
function makeDataToScreen(points: Point[], ring: Point[]) {
    const fitAxis = (axis: 0 | 1) => {
        const [i, j] = [0, 2];
        const scale = (points[j][axis] - points[i][axis]) / (ring[j][axis] - ring[i][axis]);
        const offset = points[i][axis] - ring[i][axis] * scale;
        return (v: number) => v * scale + offset;
    };
    const fx = fitAxis(0);
    const fy = fitAxis(1);
    // sanity check against the remaining vertices
    for (let k = 0; k < points.length; k++) {
        expect(Math.abs(fx(ring[k][0]) - points[k][0])).toBeLessThan(3);
        expect(Math.abs(fy(ring[k][1]) - points[k][1])).toBeLessThan(3);
    }
    return (p: Point): Point => [fx(p[0]), fy(p[1])];
}

test.describe("Scatter selection tools", () => {
    test.setTimeout(300_000);

    test("polygon, freehand, modify and pan work in place, popped out and popped back in", async ({ page }) => {
        const title = `Selection Tools ${Date.now()}`;
        let handle: SyntheticAnndataTemporaryProjectHandle | undefined;
        try {
            handle = await createTemporaryProjectViaSyntheticAnndata(page, {
                nameSegment: `scatter-selection-tools--${Date.now()}`,
                synthetic: {
                    // realistic rather than minimal for an X_umap pair to plot; kept small for software picking
                    profile: "realistic",
                    nCells: 200,
                    nGenes: 12,
                    force: true,
                },
            });

            const id = await addScatterChart(page, title);
            let box = await measureCanvas(page, id);
            await warmUp(page, page, id, box);

            await checkAllTools(page, page, id, box, "in place");

            const [popup] = await Promise.all([
                page.waitForEvent("popup"),
                page.evaluate((id) => {
                    const cm = (window as any).mdv.chartManager;
                    cm._popOutChart(cm.charts[id].chart);
                }, id),
            ]);
            try {
                await popup.waitForLoadState("domcontentloaded");
                await expect(chartScope(popup, id).locator('button[aria-label="Polygon"]')).toBeVisible();
                const popupBox = await measureCanvas(popup, id);

                await checkAllTools(popup, page, id, popupBox, "(g) popped out");
            } finally {
                if (!popup.isClosed()) {
                    await popup.close({ runBeforeUnload: true });
                }
            }

            await expect(chartScope(page, id).locator('button[aria-label="Polygon"]')).toBeVisible();
            await expect.poll(async () => (await measureCanvas(page, id)).width).toBeGreaterThan(100);
            box = await measureCanvas(page, id);

            await checkAllTools(page, page, id, box, "(h) popped back in");
        } finally {
            if (handle !== undefined) {
                await handle.cleanup();
            }
        }
    });
});
