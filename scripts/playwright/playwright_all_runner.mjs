import { spawn } from "node:child_process";
import path from "node:path";
import { fileURLToPath } from "node:url";

const HERE = path.dirname(fileURLToPath(import.meta.url));
const REPO_ROOT = path.resolve(HERE, "../..");
const args = process.argv.slice(2);

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

await run("node", ["scripts/playwright/playwright_catalog_runner.mjs", ...args]);
await run("node", ["scripts/playwright/playwright_project_runner.mjs", ...args]);
