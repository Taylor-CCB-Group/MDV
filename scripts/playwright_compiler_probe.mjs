// Quick ad-hoc probe to compare project page load between two URLs.
// Usage: node scripts/playwright_compiler_probe.mjs <url1> <url2>
// Example:
//   node scripts/playwright_compiler_probe.mjs \
//     http://127.0.0.1:5055/project/191 \
//     http://127.0.0.1:5170/project/191

import { chromium } from "@playwright/test";

const targets = process.argv.slice(2);
if (targets.length === 0) {
    console.error("Usage: node playwright_compiler_probe.mjs <url> [<url> ...]");
    process.exit(2);
}

async function probe(url) {
    let browser;
    let ctx;
    let page;
    let probeError;
    try {
        browser = await chromium.launch({ headless: true });
        ctx = await browser.newContext({ viewport: { width: 1600, height: 1000 } });
        page = await ctx.newPage();
        const issues = [];
        page.on("pageerror", (e) => issues.push({ kind: "pageerror", msg: e.message }));
        page.on("console", (m) => {
            if (m.type() === "error" || m.type() === "warning") {
                issues.push({ kind: m.type(), msg: m.text().slice(0, 400) });
            }
        });

        const t0 = Date.now();
        await page.goto(url, { waitUntil: "domcontentloaded", timeout: 60_000 });
        const tDom = Date.now();
        await page.waitForFunction(
            () => Boolean(window.mdv?.chartManager?.viewManager?.current_view),
            null,
            { timeout: 120_000 },
        );
        const tCm = Date.now();
        await page.waitForFunction(
            () => Object.keys(window.mdv?.chartManager?.charts ?? {}).length > 0,
            null,
            { timeout: 120_000 },
        );
        const tCharts = Date.now();
        await page.evaluate(() =>
            new Promise((r) => {
                let n = 4;
                const step = () => {
                    n -= 1;
                    if (n <= 0) r();
                    else requestAnimationFrame(step);
                };
                requestAnimationFrame(step);
            }),
        );
        const tSettled = Date.now();

        const summary = await page.evaluate(() => {
            const cm = window.mdv.chartManager;
            const charts = Object.values(cm.charts ?? {});
            const ds = charts.map((c) => c.chart?.dataStore).find(Boolean);
            return {
                currentView: cm.viewManager?.current_view ?? null,
                chartCount: charts.length,
                chartTypes: [...new Set(charts.map((c) => c.chart?.config?.type))],
                rowCount: ds?.size ?? null,
            };
        });
        return {
            url,
            domContentLoadedMs: tDom - t0,
            chartManagerReadyMs: tCm - t0,
            firstChartsMs: tCharts - t0,
            chartsSettledMs: tSettled - t0,
            summary,
            issues,
        };
    } catch (error) {
        probeError = error;
        throw error;
    } finally {
        if (page) {
            try {
                await page.close();
            } catch {
                // Ignore cleanup failures; preserve the original probe error if any.
            }
        }
        if (ctx) {
            try {
                await ctx.close();
            } catch {
                // Ignore cleanup failures; preserve the original probe error if any.
            }
        }
        if (browser) {
            try {
                await browser.close();
            } catch {
                // Ignore cleanup failures; preserve the original probe error if any.
            }
        }
    }
}

const results = [];
for (const url of targets) {
    try {
        results.push(await probe(url));
    } catch (e) {
        results.push({ url, error: e.message });
    }
}

console.log(JSON.stringify(results, null, 2));
