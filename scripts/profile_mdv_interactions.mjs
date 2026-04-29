import { chromium } from "@playwright/test";
import fs from "node:fs";
import path from "node:path";

function readRequiredValue(argv, index, flag) {
    const next = argv[index + 1];
    if (!next || next.startsWith("--")) {
        throw new Error(`${flag} requires a value`);
    }
    return next;
}

function readRequiredList(argv, index, flag) {
    const value = readRequiredValue(argv, index, flag);
    const values = value.split(",").map((item) => item.trim()).filter(Boolean);
    if (values.length === 0) {
        throw new Error(`${flag} requires at least one non-empty value`);
    }
    return values;
}

function readRequiredInteger(argv, index, flag) {
    const value = readRequiredValue(argv, index, flag);
    const parsed = Number(value);
    if (!Number.isInteger(parsed)) {
        throw new Error(`${flag} requires a numeric value`);
    }
    return parsed;
}

function readArgs(argv) {
    const args = {
        categories: ["type_0", "type_1", "type_2"],
        headed: false,
        out: "output/playwright/mdv-interaction-profile.json",
        repeats: 1,
        settleMs: 1000,
        url: "http://localhost:5055/project/191",
        views: ["react view", "classic view"],
    };
    for (let index = 2; index < argv.length; index += 1) {
        const arg = argv[index];
        if (arg === "--headed") {
            args.headed = true;
            continue;
        } else if (arg === "--url") {
            args.url = readRequiredValue(argv, index, arg);
        } else if (arg === "--views") {
            args.views = readRequiredList(argv, index, arg);
        } else if (arg === "--categories") {
            args.categories = readRequiredList(argv, index, arg);
        } else if (arg === "--out") {
            args.out = readRequiredValue(argv, index, arg);
        } else if (arg === "--repeats") {
            args.repeats = readRequiredInteger(argv, index, arg);
        } else if (arg === "--settle-ms") {
            args.settleMs = readRequiredInteger(argv, index, arg);
        } else {
            throw new Error(`Unknown argument: ${arg}`);
        }
        index += 1;
    }
    if (args.repeats <= 0) {
        throw new Error("--repeats must be greater than zero");
    }
    if (args.settleMs < 0) {
        throw new Error("--settle-ms must be zero or greater");
    }
    return args;
}

function viewUrl(baseUrl, view) {
    const url = new URL(baseUrl);
    url.searchParams.set("view", view);
    return url.toString();
}

async function profileView(page, baseUrl, viewName, categories, settleMs) {
    const url = viewUrl(baseUrl, viewName);
    const startedAt = Date.now();
    await page.goto(url, { waitUntil: "domcontentloaded", timeout: 60_000 });
    await page.waitForFunction(
        () => Boolean(window.mdv?.chartManager?.viewManager?.current_view),
        null,
        { timeout: 120_000 },
    );
    await page.waitForSelector(".row-bar", { timeout: 120_000 });
    await page.waitForTimeout(settleMs);
    const readyAt = Date.now();

    const metadata = await page.evaluate(() => {
        const cm = window.mdv.chartManager;
        const charts = Object.values(cm.charts).map((entry) => ({
            id: entry.chart.config.id,
            type: entry.chart.config.type,
            title: entry.chart.config.title,
            param: entry.chart.config.param,
        }));
        const rowChart = Object.values(cm.charts)
            .map((entry) => entry.chart)
            .find((chart) => chart.config.type === "row_chart" && !chart.config.wordcloud);
        const rowParam = rowChart?.config?.param?.[0] ?? null;
        const rowValues = rowParam ? rowChart.dataStore.getColumnValues(rowParam) : [];
        const dataStore = rowChart?.dataStore;
        return {
            currentView: cm.viewManager.current_view,
            charts,
            rowParam,
            rowValues,
            rowCount: dataStore?.size ?? null,
            filterSize: dataStore?.filterSize ?? null,
        };
    });

    const interactions = [];
    for (const category of categories) {
        const result = await page.evaluate(async (cat) => {
            const waitForFrames = (count) =>
                new Promise((resolve) => {
                    const step = () => {
                        count -= 1;
                        if (count <= 0) {
                            resolve();
                            return;
                        }
                        requestAnimationFrame(step);
                    };
                    requestAnimationFrame(step);
                });
            const cm = window.mdv.chartManager;
            const rowChart = Object.values(cm.charts)
                .map((entry) => entry.chart)
                .find((chart) => chart.config.type === "row_chart" && !chart.config.wordcloud);
            if (!rowChart) throw new Error("No row chart found");
            const dataStore = rowChart.dataStore;
            const label = [...document.querySelectorAll(".row-text")]
                .find((element) => element.textContent === cat);
            if (!label) throw new Error(`No row chart label found for ${cat}`);

            const before = performance.now();
            const eventPromise = new Promise((resolve, reject) => {
                const listenerId = `profile-${Date.now()}-${Math.random()}`;
                const timeoutId = setTimeout(() => {
                    dataStore.removeListener(listenerId);
                    reject(new Error(`Timed out waiting for filtered event for ${cat}`));
                }, 15_000);
                dataStore.addListener(listenerId, (type, data) => {
                    if (type !== "filtered") return;
                    dataStore.removeListener(listenerId);
                    clearTimeout(timeoutId);
                    resolve({
                        eventAt: performance.now(),
                        filterSize: dataStore.filterSize,
                        eventData: typeof data === "string" ? data : data?.constructor?.name ?? null,
                    });
                });
            });
            label.dispatchEvent(
                new MouseEvent("click", {
                    bubbles: true,
                    cancelable: true,
                    view: window,
                }),
            );
            const event = await eventPromise;
            const filteredIndicesStartedAt = performance.now();
            const filteredIndices = await dataStore.getFilteredIndices();
            const filteredIndicesAt = performance.now();
            await waitForFrames(2);
            const settledAt = performance.now();
            return {
                category: cat,
                filterEventMs: event.eventAt - before,
                filteredIndicesMs: filteredIndicesAt - filteredIndicesStartedAt,
                visualSettleMs: settledAt - before,
                filterSize: event.filterSize,
                filteredIndicesLength: filteredIndices.length,
                eventData: event.eventData,
            };
        }, category);
        interactions.push(result);
    }

    await page.evaluate(() => {
        const rowChart = Object.values(window.mdv.chartManager.charts)
            .map((entry) => entry.chart)
            .find((chart) => chart.config.type === "row_chart" && !chart.config.wordcloud);
        rowChart?.removeFilter();
    });
    await page.waitForTimeout(250);

    return {
        view: viewName,
        url,
        loadMs: readyAt - startedAt,
        metadata,
        interactions,
    };
}

const args = readArgs(process.argv);
const browser = await chromium.launch({ headless: !args.headed });
const page = await browser.newPage({ viewport: { width: 1800, height: 1100 } });
const consoleMessages = [];
page.on("console", (message) => {
    const text = message.text();
    if (/error|warn|view_loaded|gridstack/i.test(text)) {
        consoleMessages.push({
            type: message.type(),
            text: text.slice(0, 500),
        });
    }
});
page.on("pageerror", (error) => {
    consoleMessages.push({ type: "pageerror", text: error.message });
});

const runs = [];
for (let repeat = 0; repeat < args.repeats; repeat += 1) {
    for (const view of args.views) {
        runs.push(await profileView(page, args.url, view, args.categories, args.settleMs));
    }
}

await browser.close();

const profile = {
    createdAt: new Date().toISOString(),
    args,
    runs,
    consoleMessages,
};

fs.mkdirSync(path.dirname(args.out), { recursive: true });
fs.writeFileSync(args.out, `${JSON.stringify(profile, null, 2)}\n`);
console.log(JSON.stringify(profile, null, 2));
