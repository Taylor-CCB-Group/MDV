// PROTOTYPE (layout lab) — throwaway. Checks the four glitches reported on the round-one lab
// (chart leaving its card, dragging to the bottom, resizing at the bottom), for each engine.
// Usage, from the repo root: node src/react/prototypes/layoutLab/checks/glitches.mjs <projectId> [inhouse,gridstack]
// LAB_URL (default http://localhost:5170). Screenshots go to artifacts/playwright-cli/layout-lab/.
import { createRequire } from "node:module";
const require = createRequire(`${process.cwd()}/package.json`);
const { chromium } = require("@playwright/test");
const dir = `${process.cwd()}/artifacts/playwright-cli/layout-lab/`;
require("node:fs").mkdirSync(dir, { recursive: true });
const project = process.argv[2];
if (!project) throw new Error("Pass a project id");
const engines = (process.argv[3] ?? "inhouse,gridstack").split(",");
const ITEM = { inhouse: ".lab-grid > .lab-card", gridstack: ".lab-gs > .grid-stack-item" };
const HANDLE_S = { inhouse: ".lab-rz-s", gridstack: ".ui-resizable-s" };

const browser = await chromium.launch({ headless: true, args: ["--use-angle=metal", "--enable-gpu", "--ignore-gpu-blocklist"] });
const out = {};
let r0;
for (const engine of engines) {
    const page = await browser.newPage({ viewport: { width: 1440, height: 900 } });
    const errs = [];
    page.on("pageerror", (e) => errs.push(e.message.slice(0, 100)));
    r0 = errs;
    await page.goto(`${process.env.LAB_URL ?? "http://localhost:5170"}/project/${project}?layoutLab=${engine}`);
    await page.waitForSelector(ITEM[engine], { timeout: 120000 });
    await page.click('[data-act="reset"]');
    await page.waitForTimeout(4000);
    const r = {};
    const item = (i) => page.locator(ITEM[engine]).nth(i);
    const chartsInCards = () =>
        page.evaluate(() =>
            [...document.querySelectorAll(".lab-slot")].filter((s) => {
                const d = s.firstElementChild;
                return d && (Math.abs(d.offsetLeft) > 1 || Math.abs(d.offsetTop) > 1);
            }).length,
        );
    const scroller = () => page.evaluate(() => {
        const s = document.querySelector(".lab-scroll");
        return { top: Math.round(s.scrollTop), height: s.scrollHeight };
    });

    // 1. drag by the title bar: the chart must stay inside its card
    const t = await item(1).locator(".ciview-chart-titlebar").first().boundingBox();
    await page.mouse.move(t.x + 40, t.y + t.height / 2);
    await page.mouse.down();
    for (let i = 1; i <= 40; i++) { await page.mouse.move(t.x + 40 + i * 6, t.y + t.height / 2 + i * 4); await page.waitForTimeout(16); }
    await page.screenshot({ path: `${dir}glitch-${engine}-middrag.png` });
    await page.mouse.up();
    await page.waitForTimeout(800);
    r.chartsOutOfTheirCard = await chartsInCards();

    // 2. hold a dragged chart at the bottom edge: the page should scroll and grow
    await page.evaluate(() => { document.querySelector(".lab-scroll").scrollTop = 0; });
    await page.waitForTimeout(300);
    const before = await scroller();
    const t0 = await item(0).locator(".ciview-chart-titlebar").first().boundingBox();
    await page.mouse.move(t0.x + 40, t0.y + t0.height / 2);
    await page.mouse.down();
    for (let i = 1; i <= 30; i++) { await page.mouse.move(t0.x + 40, t0.y + t0.height / 2 + ((890 - t0.y) * i) / 30); await page.waitForTimeout(16); }
    for (let i = 0; i < 90; i++) { await page.mouse.move(t0.x + 40 + (i % 2), 890); await page.waitForTimeout(16); }
    const during = await scroller();
    await page.mouse.up();
    await page.waitForTimeout(800);
    r.holdAtBottom = { scrolledBy: during.top - before.top, grewBy: during.height - before.height };
    r.chartsOutOfTheirCardAfterHold = await chartsInCards();

    // 4. resize from the bottom edge, slowly, past the bottom of the window: heights should only grow
    {
        const card = item(2);
        // put the card's bottom edge 150 px above the bottom of the window
        await page.evaluate(({ sel }) => {
            const el = document.querySelectorAll(sel)[2];
            const s = document.querySelector(".lab-scroll");
            const r = el.getBoundingClientRect();
            s.scrollTop += r.bottom - (window.innerHeight - 150);
        }, { sel: ITEM[engine] });
        await page.waitForTimeout(500);
        const cb = await card.boundingBox();
        await page.mouse.move(cb.x + cb.width / 2, cb.y + cb.height / 2);
        await page.waitForTimeout(200);
        const hb = await card.locator(HANDLE_S[engine]).first().boundingBox();
        if (!hb) {
            r.resizeAtBottom = "no bottom handle";
        } else {
            const sx = hb.x + hb.width / 2, sy = hb.y + hb.height / 2;
            r.resizeStartHit = await page.evaluate(({ x, y }) => { const e = document.elementFromPoint(x, y); return e ? `${e.tagName}.${typeof e.className === "string" ? e.className : ""}` : null; }, { x: sx, y: sy });
            await page.mouse.move(sx, sy);
            await page.mouse.down();
            const heights = [];
            for (let i = 1; i <= 40; i++) { await page.mouse.move(sx, sy + (140 * i) / 40); await page.waitForTimeout(20); heights.push(Math.round((await card.boundingBox()).height)); }
            // hold at the bottom edge, jiggling by a pixel, as a hand does
            for (let i = 0; i < 60; i++) { await page.mouse.move(sx + (i % 2), 893 + (i % 3)); await page.waitForTimeout(20); heights.push(Math.round((await card.boundingBox()).height)); }
            await page.mouse.up();
            let reversals = 0;
            for (let i = 2; i < heights.length; i++) {
                const a = heights[i - 1] - heights[i - 2], b = heights[i] - heights[i - 1];
                if (a * b < 0) reversals++;
            }
            r.resizeAtBottom = { start: cb.height | 0, afterDrag: heights[39], afterHold: heights.at(-1), reversals };
        }
    }
    r.pageErrors = errs;
    out[engine] = r;
    await page.close();
}
console.log(JSON.stringify(out, null, 1));
await browser.close();
