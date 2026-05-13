import { spawn } from "node:child_process";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { hasPositionalArgs } from "./arg_utils.mjs";

const HERE = path.dirname(fileURLToPath(import.meta.url));
const REPO_ROOT = path.resolve(HERE, "../..");
const args = process.argv.slice(2);
const DEFAULT_DEV_CATALOG_SPECS = [
    "tests_playwright/catalog/create_project.spec.ts",
    "tests_playwright/catalog/import_project.spec.ts",
];

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
    forwardedArgs.unshift(...DEFAULT_DEV_CATALOG_SPECS);
}

await run("node", ["scripts/playwright/run_playwright_cli.mjs", "test", ...forwardedArgs, "--project=chromium", "--ui"]);
