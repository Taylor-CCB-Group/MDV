/**
 * End-to-end guard for the deck scatter selection tools (issues #577, #578):
 * - a double-click finishes a Polygon with the double-clicked point as its last vertex, and doesn't zoom;
 * - Freehand, Modify (vertex removal) and Pan (translate selection / controller pan + zoom) work;
 * - all of that still works after the chart is popped out into another window, and after it comes back;
 * - with nothing chosen in this browser, a plain drag moves the view and asks what dragging should
 *   do; once "Draw a selection" is chosen, a plain drag draws a rectangle, while ⌘-drag (Ctrl-drag
 *   off a Mac) and a right-button drag move the view (with a hint until the view has been moved
 *   that way), and a click on empty space clears the selection.
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

/** The selection toolbar, which starts collapsed to one button that switches to the other tool. */
function toolbar(ui: Page, id: string): Locator {
    return chartScope(ui, id).getByRole("group", { name: "choose tool for manipulating view or selection" });
}

async function selectTool(ui: Page, id: string, tool: Tool) {
    const showTools = toolbar(ui, id).locator('button[aria-label="Show tools"]');
    if (await showTools.isVisible()) await showTools.click();
    await toolbar(ui, id).locator(`button[aria-label="${tool}"]`).click();
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

/**
 * What a plain drag does in charts, for this test's browser context (`chartDragSetting`), set before
 * the app loads so that its prompt doesn't come up in tests about something else.
 */
async function chooseChartDrag(page: Page, value: "view" | "select") {
    await page.addInitScript((value) => localStorage.setItem("mdv.settings.chartDrag", value), value);
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

/** (f): with no selection, Pan (the lock) drags the view, and the wheel zooms with the view key held. */
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
    // a plain wheel is left to scroll the page; Ctrl / ⌘ + wheel zooms
    const viewKey = await viewKeyOf(ui);
    await ui.keyboard.down(viewKey);
    await ui.mouse.wheel(0, -400);
    await ui.keyboard.up(viewKey);
    await expect.poll(async () => (await readState(page, id)).viewState?.zoom).not.toBe(before.zoom);
}

/** The key the app reads to move the view: ⌘ on a Mac (where Ctrl-click is a right-click), Ctrl elsewhere. */
async function viewKeyOf(ui: Page) {
    const mac = await ui.evaluate(() => /Mac|iPhone|iPad/.test(navigator.userAgent));
    return mac ? "Meta" : "Control";
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
 * Starts recording what the chart's gesture hint shows. The hint is only up for about a second,
 * which a slow first drag (software rendering) can outlast, so it is watched rather than polled.
 */
async function watchHint(page: Page, id: string) {
    await page.evaluate((id) => {
        const content: HTMLElement = (window as any).mdv.chartManager.charts[id].chart.contentDiv;
        const hint = content.querySelector<HTMLElement>("[data-view-gesture-hint]");
        if (!hint) throw new Error("no gesture hint in the chart");
        const shown: string[] = [];
        (window as any).__hintShown = shown;
        (window as any).__hintObserver?.disconnect();
        const observer = new MutationObserver(() => {
            if (hint.style.opacity === "1") shown.push(hint.textContent ?? "");
        });
        observer.observe(hint, { attributes: true, attributeFilter: ["style"], subtree: true, characterData: true });
        (window as any).__hintObserver = observer;
    }, id);
}

/** What the hint has shown since `watchHint`. */
async function hintShown(page: Page): Promise<string[]> {
    return await page.evaluate(() => [...((window as any).__hintShown ?? [])]);
}

/** A drag from one point to another, with a button and keys held (as a real mouse and keyboard). */
async function dragWith(ui: Page, from: Point, to: Point, { button = "left", keys = [] as string[] } = {}) {
    for (const key of keys) await ui.keyboard.down(key);
    await ui.mouse.move(...from, { steps: 3 });
    await ui.mouse.down({ button: button as "left" | "right" });
    await ui.mouse.move(...to, { steps: 10 });
    await ui.mouse.up({ button: button as "left" | "right" });
    for (const key of keys) await ui.keyboard.up(key);
    await settle(ui);
}

/** Let deck and React catch up: two animation frames in the chart's window, then a short wait. */
async function settle(ui: Page) {
    await ui.evaluate(() => new Promise((r) => requestAnimationFrame(() => requestAnimationFrame(r))));
    await ui.waitForTimeout(300);
}

/** Pointer moves with no button held, straight after finishing a shape - what a real mouse does next. */
async function hoverAfterFinishing(ui: Page, box: Box) {
    for (const [fx, fy] of [
        [0.5, 0.5],
        [0.52, 0.48],
        [0.55, 0.5],
    ]) {
        await ui.mouse.move(...at(box, fx, fy));
    }
}

type SyncPointerEvent = { type: "down" | "up" | "move"; at: Point };

/**
 * Dispatch pointer events on the deck canvas in a single task, as a fast real mouse delivers them.
 * Each `page.mouse` call is its own round trip, which gives React time to re-render in between and
 * hides races between an edit and the layer getting its new data. Hover moves need the mouse event
 * as well as the pointer one: mjolnir takes hover from `mousemove`, drags from pointer events.
 */
async function dispatchInOneTask(ui: Page, id: string, events: SyncPointerEvent[]) {
    await chartScope(ui, id)
        .locator("canvas")
        .evaluateAll((canvases, events) => {
            const canvas = canvases.reduce((a, b) =>
                b.clientWidth * b.clientHeight > a.clientWidth * a.clientHeight ? b : a,
            );
            for (const { type, at } of events) {
                const buttons = type === "down" ? 1 : 0;
                const init = { clientX: at[0], clientY: at[1], button: 0, buttons, bubbles: true, cancelable: true };
                canvas.dispatchEvent(
                    new PointerEvent(`pointer${type}`, {
                        ...init,
                        pointerId: 1,
                        pointerType: "mouse",
                        isPrimary: true,
                    }),
                );
                canvas.dispatchEvent(new MouseEvent(`mouse${type}`, init));
            }
        }, events);
}

/**
 * Freehand around a square without clearing first or waiting for the commit. With `quickHover`, the
 * release and the first hover move after it arrive together, like a quick real mouse.
 */
async function dragFreehand(
    ui: Page,
    id: string,
    box: Box,
    [x0, y0, x1, y1]: [number, number, number, number],
    { quickHover }: { quickHover: boolean },
) {
    const corners = [at(box, x0, y0), at(box, x1, y0), at(box, x1, y1), at(box, x0, y1)];
    await ui.mouse.move(...corners[0], { steps: 3 });
    await ui.mouse.down();
    for (const corner of [...corners.slice(1), corners[0]]) {
        await ui.mouse.move(...corner, { steps: 8 });
    }
    if (!quickHover) {
        await ui.mouse.up();
        return;
    }
    await dispatchInOneTask(ui, id, [
        { type: "up", at: corners[0] },
        { type: "move", at: at(box, 0.5, 0.5) },
    ]);
    // keep the browser's own button state in step; the drag has already ended
    await ui.mouse.up();
}

function centreX(ring: Point[]) {
    const xs = ring.map((p) => p[0]);
    return (Math.min(...xs) + Math.max(...xs)) / 2;
}

/** The single selected ring once things have settled. */
async function settledRing(ui: Page, page: Page, id: string): Promise<Point[]> {
    await settle(ui);
    const { rings } = await readState(page, id);
    expect(rings).toHaveLength(1);
    return rings[0];
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
        await chooseChartDrag(page, "view");
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
                await expect(toolbar(popup, id)).toBeVisible();
                const popupBox = await measureCanvas(popup, id);

                await checkAllTools(popup, page, id, popupBox, "(g) popped out");
            } finally {
                if (!popup.isClosed()) {
                    await popup.close({ runBeforeUnload: true });
                }
            }

            await expect(toolbar(page, id)).toBeVisible();
            await expect.poll(async () => (await measureCanvas(page, id)).width).toBeGreaterThan(100);
            box = await measureCanvas(page, id);

            await checkAllTools(page, page, id, box, "(h) popped back in");
        } finally {
            if (handle !== undefined) {
                await handle.cleanup();
            }
        }
    });

    test("dragging moves the view until drawing is chosen; then ⌘-drag or a right-button drag moves it", async ({
        page,
    }) => {
        const title = `Drag To Select ${Date.now()}`;
        let handle: SyntheticAnndataTemporaryProjectHandle | undefined;
        try {
            handle = await createTemporaryProjectViaSyntheticAnndata(page, {
                nameSegment: `scatter-drag-to-select--${Date.now()}`,
                synthetic: { profile: "realistic", nCells: 200, nGenes: 12, force: true },
            });
            const id = await addScatterChart(page, title);
            const box = await measureCanvas(page, id);
            const viewKey = await viewKeyOf(page);
            const mac = viewKey === "Meta";
            // a fresh browser context: nothing chosen, and no view moved with a modifier yet
            let start = await waitForViewState(page, id);

            await test.step("with nothing chosen, a plain drag moves the view and asks what dragging should do", async () => {
                await dragWith(page, at(box, 0.3, 0.3), at(box, 0.45, 0.4));
                await expect.poll(async () => (await readState(page, id)).viewState?.target).not.toEqual(start.target);
                expect((await readState(page, id)).rings).toHaveLength(0);
                // collapsed, the toolbar offers drawing as the other tool
                await expect(toolbar(page, id).locator('button[aria-label="Rectangle"]')).toBeVisible();
                const choice = chartScope(page, id).getByRole("dialog", { name: "What dragging does on charts" });
                await expect(choice).toBeVisible();
                await choice.getByRole("button", { name: "Draw a selection" }).click();
                await expect(choice).toBeHidden();
                expect(await page.evaluate(() => localStorage.getItem("mdv.settings.chartDrag"))).toBe("select");
                // the chart now draws, and its toolbar offers moving the view instead
                await expect(toolbar(page, id).locator('button[aria-label="Pan"]')).toBeVisible();
                start = await waitForViewState(page, id);
            });

            await test.step("then a plain drag draws a rectangle, and says how to move the view", async () => {
                await watchHint(page, id);
                await page.mouse.move(...at(box, 0.3, 0.3), { steps: 3 });
                await page.mouse.down();
                await page.mouse.move(...at(box, 0.6, 0.6), { steps: 10 });
                await page.mouse.up();
                const ring = await settledRing(page, page, id);
                await expect
                    .poll(() => hintShown(page))
                    .toContainEqual(`Hold ${mac ? "⌘" : "Ctrl"} and drag to move the view`);
                expect(ring).toHaveLength(5);
                expect((await readState(page, id)).viewState?.target).toEqual(start.target);
            });

            const drawn = (await readState(page, id)).rings[0];
            let view = await waitForViewState(page, id);

            await test.step(`${viewKey}-drag moves the view and leaves the selection`, async () => {
                await dragWith(page, at(box, 0.7, 0.2), at(box, 0.8, 0.35), { keys: [viewKey] });
                const { rings, viewState } = await readState(page, id);
                expect(viewState?.target).not.toEqual(view.target);
                expect(rings).toEqual([drawn]);
                expect(await page.evaluate(() => localStorage.getItem("mdv.viewDragLearned"))).toBe("1");
                view = await waitForViewState(page, id);
            });

            await test.step("a right-button drag moves the view too", async () => {
                await dragWith(page, at(box, 0.7, 0.2), at(box, 0.6, 0.3), { button: "right" });
                const { rings, viewState } = await readState(page, id);
                expect(viewState?.target).not.toEqual(view.target);
                expect(rings).toEqual([drawn]);
                view = await waitForViewState(page, id);
            });

            await test.step("Shift-drag still draws, and the hint has stopped", async () => {
                await watchHint(page, id);
                await page.keyboard.down("Shift");
                await page.mouse.move(...at(box, 0.2, 0.6), { steps: 3 });
                await page.mouse.down();
                await page.mouse.move(...at(box, 0.4, 0.8), { steps: 10 });
                await page.mouse.up();
                await page.keyboard.up("Shift");
                const ring = await settledRing(page, page, id);
                expect(await hintShown(page)).toEqual([]);
                expect(ring).not.toEqual(drawn);
                expect((await readState(page, id)).viewState?.target).toEqual(view.target);
            });

            await test.step("a click on empty space clears the selection", async () => {
                // zoomed well out, the data and the selection shrink to the middle and leave the edges empty
                await page.evaluate((id) => {
                    const config = (window as any).mdv.chartManager.charts[id].chart.config;
                    const vs = config.viewState;
                    const out = (z: number | undefined) => (z === undefined ? undefined : z - 3);
                    config.viewState = { ...vs, zoom: out(vs.zoom), zoomX: out(vs.zoomX), zoomY: out(vs.zoomY) };
                }, id);
                await settle(page);
                expect((await readState(page, id)).rings).toHaveLength(1);
                // the left edge, half way down (the bottom of a 700px chart is out of the window)
                await page.mouse.click(...at(box, 0.02, 0.5));
                await expect.poll(async () => (await readState(page, id)).rings.length).toBe(0);
            });
        } finally {
            if (handle !== undefined) {
                await handle.cleanup();
            }
        }
    });

    /**
     * The layer only gets a new selection when React re-renders, and until then editable-layers hands
     * edits the previous one: writing those back put the old shape back as soon as the pointer moved.
     */
    test("a shape drawn over an existing selection stays selected when the pointer moves on", async ({ page }) => {
        const title = `Selection Persistence ${Date.now()}`;
        let handle: SyntheticAnndataTemporaryProjectHandle | undefined;
        await chooseChartDrag(page, "view");
        try {
            handle = await createTemporaryProjectViaSyntheticAnndata(page, {
                nameSegment: `scatter-selection-persistence--${Date.now()}`,
                synthetic: { profile: "realistic", nCells: 200, nGenes: 12, force: true },
            });
            const id = await addScatterChart(page, title);
            const box = await measureCanvas(page, id);
            await warmUp(page, page, id, box);

            await selectTool(page, id, "Freehand");
            await dragFreehand(page, id, box, [0.15, 0.3, 0.4, 0.7], { quickHover: false });
            await expect.poll(async () => (await readState(page, id)).rings.length).toBe(1);
            let previous = await settledRing(page, page, id);

            await test.step("freehand over freehand", async () => {
                await dragFreehand(page, id, box, [0.6, 0.3, 0.85, 0.7], { quickHover: true });
                await hoverAfterFinishing(page, box);
                const ring = await settledRing(page, page, id);
                expect(ring).not.toEqual(previous);
                expect(centreX(ring)).toBeGreaterThan(centreX(previous));
                previous = ring;
            });

            await test.step("polygon over freehand", async () => {
                await selectTool(page, id, "Polygon");
                const points = [at(box, 0.15, 0.3), at(box, 0.4, 0.3), at(box, 0.4, 0.7), at(box, 0.2, 0.75)];
                for (const [x, y] of points.slice(0, 3)) {
                    await page.mouse.move(x, y, { steps: 3 });
                    await page.mouse.click(x, y);
                    await pause(page);
                }
                await page.mouse.move(...points[3], { steps: 3 });
                // double-click the last vertex and move on, in one task
                await dispatchInOneTask(page, id, [
                    { type: "down", at: points[3] },
                    { type: "up", at: points[3] },
                    { type: "down", at: points[3] },
                    { type: "up", at: points[3] },
                    { type: "move", at: at(box, 0.5, 0.5) },
                ]);
                await hoverAfterFinishing(page, box);
                const ring = await settledRing(page, page, id);
                expect(ring).toHaveLength(5);
                expect(centreX(ring)).toBeLessThan(centreX(previous));
                previous = ring;
            });

            await test.step("quick rectangle over polygon", async () => {
                await selectTool(page, id, "Rectangle");
                await page.mouse.move(...at(box, 0.6, 0.3), { steps: 3 });
                await page.mouse.down();
                await page.mouse.move(...at(box, 0.85, 0.7));
                await page.mouse.up();
                await hoverAfterFinishing(page, box);
                const ring = await settledRing(page, page, id);
                expect(ring).toHaveLength(5);
                expect(centreX(ring)).toBeGreaterThan(centreX(previous));
            });
        } finally {
            if (handle !== undefined) {
                await handle.cleanup();
            }
        }
    });
});
