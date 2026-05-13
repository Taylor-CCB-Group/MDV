import fs from "node:fs/promises";
import os from "node:os";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { execFile } from "node:child_process";
import { promisify } from "node:util";

const execFileAsync = promisify(execFile);
const HERE = path.dirname(fileURLToPath(import.meta.url));
const REPO_ROOT = path.resolve(HERE, "../..");
const TEST_BASE_URL = process.env.TEST_BASE_URL || "http://localhost:5055";
const DIAGNOSTIC = process.argv.includes("--diagnostic");
let pythonBin;

function fail(message) {
    console.error(`[playwright project preflight] ${message}`);
    process.exit(1);
}

async function ensurePlaywrightInstalled() {
    try {
        await import("@playwright/test");
    } catch (error) {
        const message = error instanceof Error ? error.message : String(error);
        fail(`@playwright/test is unavailable. Run \`pnpm install\`. ${message}`);
    }
}

async function resolvePythonVenv() {
    const candidates = [
        path.join(REPO_ROOT, "venv", "bin", "python"),
        path.join(REPO_ROOT, "python", ".venv", "bin", "python"),
    ];

    for (const candidate of candidates) {
        try {
            await fs.access(candidate);
            return candidate;
        } catch {
            // try next
        }
    }

    fail(
        `Python venv is missing. Looked for ${candidates.join(" and ")}. Recreate it with \`pnpm run python-setup\` or your documented local Python setup.`,
    );
}

async function ensureMdvRoot() {
    const mdvRoot = path.join(os.homedir(), "mdv");
    await fs.mkdir(mdvRoot, { recursive: true });
}

async function ensureBackendReachable() {
    let response;
    try {
        response = await fetch(TEST_BASE_URL, { method: "GET" });
    } catch (error) {
        const message = error instanceof Error ? error.message : String(error);
        fail(`Backend is unreachable at ${TEST_BASE_URL}. ${message}`);
    }

    if (!response.ok) {
        fail(`Backend check failed at ${TEST_BASE_URL} with status ${response.status}.`);
    }
}

async function runDiagnosticPythonImports() {
    const command = [
        "-c",
        [
            "import h5py",
            "import numpy",
            "import polars",
            "import scanpy",
            "print('ok')",
        ].join("\n"),
    ];

    try {
        await execFileAsync(pythonBin, command, {
            cwd: path.join(REPO_ROOT, "python"),
            env: {
                ...process.env,
                PYTHONPATH: path.join(REPO_ROOT, "python"),
            },
        });
    } catch (error) {
        const message = error instanceof Error ? error.message : String(error);
        fail(
            `Diagnostic Python import check failed for the synthetic AnnData fixture path. ${message}`,
        );
    }
}

await ensurePlaywrightInstalled();
pythonBin = await resolvePythonVenv();
await ensureMdvRoot();
await ensureBackendReachable();

if (DIAGNOSTIC) {
    await runDiagnosticPythonImports();
}

console.log(
    `[playwright project preflight] ok (${DIAGNOSTIC ? "diagnostic" : "light"})`,
);
