import { spawnSync } from "node:child_process";
import { existsSync } from "node:fs";
import path from "node:path";
import { pathToFileURL } from "node:url";

const outputPath = process.env.VITE_BUNDLE_ANALYZE_OUTPUT || "dist/bundle-analysis.html";
const absolutePath = path.resolve(process.cwd(), outputPath);

if (!existsSync(absolutePath)) {
    console.warn(`Bundle report not found at: ${absolutePath}`);
    process.exit(0);
}

const reportUrl = pathToFileURL(absolutePath).href;
const isMac = process.platform === "darwin";
const isWindows = process.platform === "win32";

const command = isMac ? "open" : isWindows ? "cmd" : "xdg-open";
const args = isMac ? [reportUrl] : isWindows ? ["/c", "start", "", reportUrl] : [reportUrl];

const result = spawnSync(command, args, { stdio: "ignore" });
if (result.status !== 0) {
    console.warn(`Could not automatically open bundle report. Open manually: ${reportUrl}`);
} else {
    console.log(`Opened bundle report: ${reportUrl}`);
}
