import { spawn } from "node:child_process";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { hasPositionalArgs } from "./arg_utils.mjs";

const HERE = path.dirname(fileURLToPath(import.meta.url));
const REPO_ROOT = path.resolve(HERE, "../..");
const args = process.argv.slice(2);
const DEFAULT_CATALOG_SPECS = [
    "tests_playwright/catalog/catalog_view.spec.ts",
    "tests_playwright/catalog/project_operations.spec.ts",
];

function hasFlag(prefix) {
    return args.some((arg) => arg === prefix || arg.startsWith(`${prefix}=`));
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
if (!hasPositionalArgs(args)) {
    forwardedArgs.unshift(...DEFAULT_CATALOG_SPECS);
}
if (!hasFlag("--project")) {
    forwardedArgs.push("--project=chromium");
}

await run("node", ["scripts/playwright/run_playwright_cli.mjs", "test", ...forwardedArgs]);
