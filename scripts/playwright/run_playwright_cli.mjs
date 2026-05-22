import fs from "node:fs/promises";
import path from "node:path";
import { spawn } from "node:child_process";
import { fileURLToPath } from "node:url";

const HERE = path.dirname(fileURLToPath(import.meta.url));
const REPO_ROOT = path.resolve(HERE, "../..");
const args = process.argv.slice(2);

function getLocalPlaywrightBinary() {
    const binaryName = process.platform === "win32" ? "playwright.cmd" : "playwright";
    return path.join(REPO_ROOT, "node_modules", ".bin", binaryName);
}

async function resolvePlaywrightCommand() {
    const localBinary = getLocalPlaywrightBinary();
    try {
        await fs.access(localBinary);
        return {
            command: localBinary,
            commandArgs: args,
        };
    } catch {
        const pnpmCommand = process.platform === "win32" ? "pnpm.cmd" : "pnpm";
        return {
            command: pnpmCommand,
            commandArgs: ["exec", "playwright", ...args],
        };
    }
}

const { command, commandArgs } = await resolvePlaywrightCommand();

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
        reject(
            new Error(
                `Playwright CLI failed with code ${code}. Run \`pnpm install\` and \`pnpm exec playwright install --with-deps\` if Playwright is missing.`,
            ),
        );
    });
    child.on("error", (error) => {
        reject(
            new Error(
                `Unable to launch Playwright CLI: ${error instanceof Error ? error.message : String(error)}`,
            ),
        );
    });
});
