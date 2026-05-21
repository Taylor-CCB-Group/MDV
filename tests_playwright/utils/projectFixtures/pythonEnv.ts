import { execFile } from "node:child_process";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { promisify } from "node:util";

const execFileAsync = promisify(execFile);

const THIS_DIR = path.dirname(fileURLToPath(import.meta.url));
/** Repository root (this file lives in `tests_playwright/utils/projectFixtures/`). */
export const REPO_ROOT = path.resolve(THIS_DIR, "..", "..", "..");

export async function getPythonCandidates(): Promise<string[]> {
    const poetryProjectDir = path.join(REPO_ROOT, "python");
    let poetryVenvPython: string | undefined;
    try {
        const { stdout } = await execFileAsync("poetry", ["env", "info", "--path"], {
            cwd: poetryProjectDir,
        });
        const poetryEnvPath = stdout.trim();
        if (poetryEnvPath) {
            poetryVenvPython = path.join(poetryEnvPath, "bin", "python");
        }
    } catch {
        // Best-effort Poetry env detection; fallback candidates below still apply.
    }

    const poetryInProjectPython = path.join(REPO_ROOT, "python", ".venv", "bin", "python");
    const repoVenvPython = path.join(REPO_ROOT, "venv", "bin", "python");
    const candidates = [
        poetryVenvPython,
        poetryInProjectPython,
        repoVenvPython,
        process.env.PYTHON,
        "python3",
        "python",
    ];
    return candidates.filter((candidate): candidate is string => Boolean(candidate));
}

export function buildPythonPath(): string {
    const repoPythonPath = path.join(REPO_ROOT, "python");
    const existing = process.env.PYTHONPATH;
    return existing ? `${repoPythonPath}${path.delimiter}${existing}` : repoPythonPath;
}
