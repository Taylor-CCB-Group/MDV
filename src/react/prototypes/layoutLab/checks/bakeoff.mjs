// PROTOTYPE (layout lab) bake-off script — throwaway, kept as a reference for the showcase's checks.
// Usage, from the repo root with the dev server and backend running:
//   LAB_PROJECTS='[{"key":"spatial","id":1},{"key":"many","id":2,"view":"Some view"}]' \
//   node src/react/prototypes/layoutLab/checks/bakeoff.mjs [smoke]
// LAB_URL (default http://localhost:5170), LAB_ENGINES (default inhouse,gridstack).
// Project ids and view names stay out of the repo: pass them in LAB_PROJECTS.
// Results and screenshots go to artifacts/playwright-cli/layout-lab/bakeoff (gitignored).
// For each project and engine: open the converted view, check the layout at four widths,
// drag and resize with real mouse input paced at ~60Hz, save, reload and compare.
import { chromium } from "@playwright/test";
import fs from "node:fs";
import path from "node:path";

const base = process.env.LAB_URL ?? "http://localhost:5170";
const out = path.resolve("artifacts/playwright-cli/layout-lab/bakeoff");
fs.mkdirSync(out, { recursive: true });
const smoke = process.argv[2] === "smoke";
const PROJECTS = JSON.parse(process.env.LAB_PROJECTS ?? "[]");
if (!PROJECTS.length) throw new Error("Set LAB_PROJECTS to a JSON list of {key, id, view?}");
const ENGINES = process.env.LAB_ENGINES?.split(",") ?? ["inhouse", "gridstack"];
const ITEM = { inhouse: ".lab-grid > .lab-card", gridstack: ".lab-gs > .grid-stack-item" };
const HANDLE_E = { inhouse: ".lab-rz-e", gridstack: ".ui-resizable-e" };
const WIDTHS = [1440, 1024, 768, 390];
const results = [];

const browser = await chromium.launch({
    headless: true,
    args: ["--use-angle=metal", "--enable-gpu", "--ignore-gpu-blocklist"],
});

async function open(page, p, engine, { clear = false } = {}) {
    await page.goto(`${base}/project/${p.id}?layoutLab=${engine}`);
    await page.waitForSelector(ITEM[engine], { timeout: 120_000 });
    if (p.view) {
        const current = await page.evaluate(() => window.mdv.chartManager.viewManager.current_view);
        if (current !== p.view) {
            await page.evaluate((v) => window.mdv.chartManager.viewManager.changeView(v), p.view);
        }
        await page.waitForFunction(
            (v) => window.mdv.chartManager.viewManager.current_view === v && !window.mdv.chartManager._toLoadCharts,
            p.view,
            { timeout: 120_000 },
        );
    }
    if (clear) {
        await page.click('[data-act="reset"]');
    }
    await page.waitForTimeout(5000); // charts finish drawing, transitions settle
}

/** Overlaps and items past the grid edge, from the DOM, plus the column count. */
async function check(page) {
    return page.evaluate(() => {
        const cards = [...document.querySelectorAll(".lab-grid > .lab-card, .lab-gs > .grid-stack-item, .lab-rgl .react-grid-item")].map((el) => {
            const r = el.getBoundingClientRect();
            return { id: el.dataset.id ?? el.getAttribute("gs-id"), l: r.left, r: r.right, t: r.top, b: r.bottom, w: r.width, h: r.height };
        });
        const grid = document.querySelector(".lab-grid, .lab-gs, .lab-rgl .react-grid-layout").getBoundingClientRect();
        const bad = [];
        for (const a of cards) {
            if (a.r > grid.right + 1 || a.l < grid.left - 1) bad.push(`${a.id} outside`);
            if (a.w < 40 || a.h < 40) bad.push(`${a.id} too small ${Math.round(a.w)}x${Math.round(a.h)}`);
            for (const b of cards)
                if (a.id < b.id && a.l < b.r - 1 && b.l < a.r - 1 && a.t < b.b - 1 && b.t < a.b - 1) bad.push(`${a.id} overlaps ${b.id}`);
        }
        // charts whose own div did not follow its slot
        for (const el of document.querySelectorAll(".lab-root .ciview-chart-panel, .lab-root [id^='chart']")) {
            void el;
        }
        const gridEl = document.querySelector(".lab-grid, .lab-gs, .lab-rgl");
        const gsCols = [...gridEl.classList].map((c) => /^gs-(\d+)$/.exec(c)?.[1]).find(Boolean);
        const cols = Number(gridEl.dataset.cols ?? gsCols ?? 12);
        const outOfCard = [...document.querySelectorAll(".lab-slot")].filter((s) => {
            const d = s.firstElementChild;
            return d && (Math.abs(d.offsetLeft) > 1 || Math.abs(d.offsetTop) > 1);
        }).length;
        if (outOfCard) bad.push(`${outOfCard} charts out of their card`);
        return { cols, n: cards.length, bad, height: Math.round(grid.height) };
    });
}

async function pacedDrag(page, from, to, steps = 90) {
    await page.mouse.move(from.x, from.y);
    await page.mouse.down();
    for (let i = 1; i <= steps; i++) {
        await page.mouse.move(from.x + ((to.x - from.x) * i) / steps, from.y + ((to.y - from.y) * i) / steps);
        await page.waitForTimeout(16);
    }
    await page.mouse.up();
}

async function lastMeasurement(page) {
    await page.waitForTimeout(900);
    return page.evaluate(() => (window.__layoutLab ?? []).at(-1));
}

const item = (page, engine, index) => page.locator(ITEM[engine]).nth(index);

for (const p of smoke ? PROJECTS.slice(0, 1) : PROJECTS) {
    for (const engine of ENGINES) {
        const page = await browser.newPage({ viewport: { width: 1440, height: 900 } });
        const errors = [];
        let phase = "load";
        page.on("pageerror", (e) => errors.push(`${phase}: ${e.message.slice(0, 120)}`));
        page.on("console", (m) => {
            if (m.type() === "error" && !m.text().includes("404")) errors.push(`${phase}: ${m.text().slice(0, 120)}`);
        });
        await page.addInitScript(() => {
            window.localStorage.setItem("MDV_SPATIAL_PERF", "1");
        });
        const r = { project: p.key, engine, widths: {}, gestures: [], errors };
        await open(page, p, engine, { clear: true });
        r.converted = await page.evaluate(() => window.__labConverted);

        for (const w of WIDTHS) {
            phase = `width ${w}`;
            await page.setViewportSize({ width: w, height: 900 });
            await page.waitForTimeout(1500);
            r.widths[w] = await check(page);
            await page.screenshot({ path: path.join(out, `${p.key}-${engine}-${w}.png`) });
        }
        await page.setViewportSize({ width: 1440, height: 900 });
        await page.waitForTimeout(1500);

        if (!smoke) {
            await page.evaluate(() => {
                window.__layoutLab = [];
                window.__resetSpatialPerf?.();
            });
            // drag the second chart by its title bar to the right and down
            phase = "move";
            await item(page, engine, 1).scrollIntoViewIfNeeded();
            await page.waitForTimeout(300);
            const t1 = await item(page, engine, 1).locator(".ciview-chart-titlebar").first().boundingBox();
            const from = { x: t1.x + 30, y: t1.y + t1.height / 2 };
            r.moveStartHit = await page.evaluate(({ x, y }) => {
                const el = document.elementFromPoint(x, y);
                return el ? `${el.tagName}.${typeof el.className === "string" ? el.className : ""}` : null;
            }, from);
            await pacedDrag(page, from, { x: from.x + 300, y: from.y + 200 });
            r.gestures.push({ kind: "move", ...(await lastMeasurement(page)) });
            await page.evaluate(() => {
                window.__layoutLab = [];
            });
            r.afterMove = await check(page);
            // resize the first chart from its right edge
            phase = "resize";
            await item(page, engine, 0).scrollIntoViewIfNeeded();
            await page.waitForTimeout(300);
            {
                const ib = await item(page, engine, 0).boundingBox();
                await page.mouse.move(ib.x + 60, Math.max(ib.y, 0) + 60);
            }
            await page.waitForTimeout(300);
            const hb = await item(page, engine, 0).locator(HANDLE_E[engine]).first().boundingBox();
            if (!hb) {
                r.gestures.push({ kind: "resize", missed: "handle not visible" });
            }
            const handle = hb
                ? { x: hb.x + hb.width / 2, y: Math.max(hb.y + 20, Math.min(hb.y + hb.height / 2, 700)) }
                : { x: 0, y: 0 };
            if (hb) {
            r.resizeHandleHit = await page.evaluate(
                ({ x, y }) => document.elementFromPoint(x, y)?.className ?? null,
                handle,
            );
            await page.mouse.move(handle.x, handle.y);
            await page.waitForTimeout(200);
            await pacedDrag(page, handle, { x: handle.x + 200, y: handle.y });
            r.gestures.push({ kind: "resize", ...((await lastMeasurement(page)) ?? { missed: true }) });
            }
            r.afterResize = await check(page);
            r.spatialPerf = await page.evaluate(() => {
                const s = window.__spatialPerf;
                if (!s) return null;
                const o = {};
                for (const [k, v] of Object.entries(s)) if (/loadRenderData|render:spatial/.test(k)) o[k] = v.count ?? v;
                return o;
            });

            // save, reload, compare
            phase = "reload";
            await page.click('[data-act="save"]');
            const before = await page.evaluate(() => window.__labLayout);
            await open(page, p, engine);
            const after = await page.evaluate(() => ({ layout: window.__labLayout, saved: window.__labSaved }));
            await page.waitForTimeout(1000);
            const settled = await page.evaluate(() => window.__labLayout);
            const same = (a, b) => JSON.stringify(Object.entries(a).sort()) === JSON.stringify(Object.entries(b).sort());
            r.reload = {
                savedMatches: same(before, after.saved ?? {}),
                loadedMatches: same(before, after.layout ?? {}),
                stillMatchesAfterSettle: same(before, settled ?? {}),
            };
            await page.screenshot({ path: path.join(out, `${p.key}-${engine}-reloaded.png`) });
            await page.click('[data-act="reset"]');
        }
        results.push(r);
        console.log(JSON.stringify({ ...r, converted: undefined }, null, 1).slice(0, 3000));
        await page.close();
    }
}
fs.writeFileSync(path.join(out, smoke ? "smoke.json" : process.env.BAKEOFF_ENGINES ? `results-${process.env.BAKEOFF_ENGINES}.json` : "results.json"), JSON.stringify(results, null, 2));
await browser.close();
