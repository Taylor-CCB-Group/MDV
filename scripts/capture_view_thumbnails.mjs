#!/usr/bin/env node
/**
 * Save each view in a running MDV project so the app stores a thumbnail.
 *
 * A thumbnail is the screenshot MDV takes in saveView (Ctrl+S) and writes to
 * views.json as viewImage. This opens the project in Chromium, switches to
 * each listed view, waits until its charts have loaded, and saves it.
 *
 * The project must already be served with edit permission. datafile.h5 is not
 * opened by this script.
 *
 *   node scripts/capture_view_thumbnails.mjs --url http://localhost:5055/project/12
 *   node scripts/capture_view_thumbnails.mjs --base-url http://localhost:5055 --project 12
 *
 * Views that already have a viewImage are skipped. --all recaptures every view.
 * --only "Cell atlas" limits the run to named views (repeat the flag to name more).
 */

import { chromium } from "playwright";

const DEFAULT_SETTLE_MS = 2000;
const DEFAULT_TIMEOUT_MS = 120000;

function usage() {
    return `Usage: node scripts/capture_view_thumbnails.mjs (--url <project-url> | --base-url <origin> --project <id>) [options]

Options:
  --only <view name>   Capture this view. Repeat for more than one.
  --all                Recapture views that already have a thumbnail.
  --settle-ms <n>      Wait after the view loads before the screenshot (default ${DEFAULT_SETTLE_MS}).
  --timeout-ms <n>     How long to wait for a view to load or save (default ${DEFAULT_TIMEOUT_MS}).
  --headed             Show the browser window.
  --help               Show this message.`;
}

function parseArgs(argv) {
    const args = {
        url: "",
        baseUrl: "",
        project: "",
        only: [],
        all: false,
        settleMs: DEFAULT_SETTLE_MS,
        timeoutMs: DEFAULT_TIMEOUT_MS,
        headed: false,
    };
    for (let i = 0; i < argv.length; i += 1) {
        const flag = argv[i];
        const value = argv[i + 1];
        if (flag === "--help") {
            console.log(usage());
            process.exit(0);
        } else if (flag === "--url") {
            args.url = value ?? "";
            i += 1;
        } else if (flag === "--base-url") {
            args.baseUrl = value ?? "";
            i += 1;
        } else if (flag === "--project") {
            args.project = value ?? "";
            i += 1;
        } else if (flag === "--only") {
            args.only.push(value ?? "");
            i += 1;
        } else if (flag === "--all") {
            args.all = true;
        } else if (flag === "--settle-ms") {
            args.settleMs = Number(value);
            i += 1;
        } else if (flag === "--timeout-ms") {
            args.timeoutMs = Number(value);
            i += 1;
        } else if (flag === "--headed") {
            args.headed = true;
        } else {
            throw new Error(`Unknown argument ${flag}\n${usage()}`);
        }
    }
    if (!args.url) {
        if (!args.baseUrl || !args.project) {
            throw new Error(`Pass --url, or both --base-url and --project.\n${usage()}`);
        }
        const origin = args.baseUrl.replace(/\/+$/, "");
        args.url = `${origin}/project/${encodeURIComponent(args.project)}`;
    }
    if (!Number.isFinite(args.settleMs) || args.settleMs < 0) {
        throw new Error("--settle-ms must be a number of milliseconds");
    }
    if (!Number.isFinite(args.timeoutMs) || args.timeoutMs <= 0) {
        throw new Error("--timeout-ms must be a positive number of milliseconds");
    }
    args.only = args.only.filter(Boolean);
    return args;
}

async function waitForApp(page, timeoutMs) {
    await page.waitForFunction(() => Boolean(window.mdv?.chartManager?.viewManager), null, { timeout: timeoutMs });
    return page.evaluate(() => {
        const config = window.mdv.chartManager.config;
        return { permission: config?.permission ?? "", views: window.mdv.chartManager.viewManager.all_views ?? [] };
    });
}

function installLoadListener(page) {
    return page.evaluate(() => {
        window.__mdvThumb = { loaded: null };
        window.mdv.chartManager.addListener("thumbnail-capture", (type, _cm, data) => {
            if (type === "view_loaded") window.__mdvThumb.loaded = data;
        });
    });
}

async function openView(page, name, timeoutMs) {
    const ready = await page.evaluate((viewName) => {
        const cm = window.mdv.chartManager;
        return cm.viewManager.current_view === viewName && !cm._toLoadCharts;
    }, name);
    if (ready) return;
    await page.evaluate((viewName) => {
        window.__mdvThumb.loaded = null;
        window.mdv.chartManager.viewManager.checkAndChangeView(viewName);
    }, name);
    await page.waitForFunction((viewName) => window.__mdvThumb.loaded === viewName, name, { timeout: timeoutMs });
}

async function viewsToCapture(page, only, captureAll) {
    const listed = await page.evaluate(async () => {
        const cm = window.mdv.chartManager;
        const names = cm.viewManager.all_views ?? [];
        const rows = [];
        for (const name of names) {
            const data = await cm.viewLoader(name);
            rows.push({ name, hasImage: Boolean(data?.viewImage) });
        }
        return rows;
    });
    const wanted = only.length > 0 ? listed.filter((row) => only.includes(row.name)) : listed;
    const missing = only.filter((name) => !listed.some((row) => row.name === name));
    for (const name of missing) {
        console.log(`skip: ${name}: not in the view list`);
    }
    return wanted.filter((row) => captureAll || !row.hasImage);
}

async function captureView(page, name, settleMs, timeoutMs) {
    await openView(page, name, timeoutMs);
    if (settleMs > 0) await page.waitForTimeout(settleMs);
    const saveResponse = page.waitForResponse(
        (response) => response.url().includes("/save_state") && response.request().method() === "POST",
        { timeout: timeoutMs },
    );
    const image = await page.evaluate(async () => {
        await window.mdv.chartManager.viewManager.saveView();
        const saved = window.mdv.chartManager.getState()?.view?.viewImage ?? "";
        return { ok: saved.startsWith("data:image"), length: saved.length };
    });
    if (!image.ok) {
        throw new Error("saveView did not produce an image");
    }
    const response = await saveResponse;
    const body = await response.json().catch(() => null);
    if (!response.ok() || body?.success === false) {
        throw new Error(`save_state failed (${response.status()})`);
    }
    return image.length;
}

async function main() {
    const args = parseArgs(process.argv.slice(2));
    const browser = await chromium.launch({ headless: !args.headed });
    const page = await browser.newPage({ viewport: { width: 1440, height: 900 } });
    try {
        console.log(`Opening ${args.url}`);
        await page.goto(args.url, { waitUntil: "domcontentloaded", timeout: args.timeoutMs });
        const app = await waitForApp(page, args.timeoutMs);
        if (app.permission !== "edit") {
            throw new Error(
                `This project is ${app.permission || "not editable"}. Thumbnails are saved only with edit permission.`,
            );
        }
        await installLoadListener(page);
        const targets = await viewsToCapture(page, args.only, args.all);
        if (targets.length === 0) {
            console.log("Nothing to capture.");
            return;
        }
        console.log(`Capturing ${targets.length} view${targets.length === 1 ? "" : "s"}`);
        let failed = 0;
        for (const target of targets) {
            try {
                const bytes = await captureView(page, target.name, args.settleMs, args.timeoutMs);
                console.log(`saved: ${target.name} (${bytes} chars)`);
            } catch (error) {
                failed += 1;
                console.log(`failed: ${target.name}: ${error instanceof Error ? error.message : error}`);
            }
        }
        if (failed > 0) process.exitCode = 1;
    } finally {
        await browser.close();
    }
}

main().catch((error) => {
    console.error(error instanceof Error ? error.message : error);
    process.exit(1);
});
