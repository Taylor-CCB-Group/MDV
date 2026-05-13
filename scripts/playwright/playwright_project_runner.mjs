import { spawn } from "node:child_process";
import path from "node:path";
import { fileURLToPath } from "node:url";

const HERE = path.dirname(fileURLToPath(import.meta.url));
const REPO_ROOT = path.resolve(HERE, "../..");
const args = process.argv.slice(2);
const FLAGS_WITH_VALUES = new Set([
    "--project",
    "-p",
    "--workers",
    "-j",
    "--grep",
    "-g",
    "--grep-invert",
    "--reporter",
    "--config",
]);

function hasFlag(prefix) {
    return args.some((arg, index) => arg === prefix || arg.startsWith(`${prefix}=`) || (arg === prefix && index < args.length - 1));
}

function hasPositionalArgs() {
    for (let index = 0; index < args.length; index += 1) {
        const arg = args[index];
        if (arg === "--") {
            return index < args.length - 1;
        }
        if (!arg.startsWith("-")) {
            return true;
        }
        if (FLAGS_WITH_VALUES.has(arg)) {
            index += 1;
        }
    }
    return false;
}

async function run(command, commandArgs) {
    await new Promise((resolve, reject) => {
        const child = spawn(command, commandArgs, {
            cwd: REPO_ROOT,
            env: process.env,
            stdio: "inherit",
        });
        child.on("exit", (code) => {
            if (code === 0) {
                resolve();
                return;
            }
            reject(new Error(`${command} exited with code ${code}`));
        });
        child.on("error", reject);
    });
}

const forwardedArgs = [...args];
if (!hasPositionalArgs()) {
    forwardedArgs.unshift("tests_playwright/project/");
}
if (!hasFlag("--project")) {
    forwardedArgs.push("--project=chromium");
}
if (!hasFlag("--workers")) {
    forwardedArgs.push("--workers=1");
}

await run("node", ["scripts/playwright/playwright_project_preflight.mjs"]);
await run("node", ["scripts/playwright/run_playwright_cli.mjs", "test", ...forwardedArgs]);
