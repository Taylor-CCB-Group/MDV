import path from "node:path";
import { fileURLToPath } from "node:url";

const THIS_DIR = path.dirname(fileURLToPath(import.meta.url));
/** Repository root (this file lives in `tests_playwright/utils/projectFixtures/`). */
export const REPO_ROOT = path.resolve(THIS_DIR, "..", "..", "..");

export async function getPythonCandidates(): Promise<string[]> {
    const pythonProjectDir = path.join(REPO_ROOT, "python");
    const uvVenvPython = path.join(pythonProjectDir, ".venv", "bin", "python");
    const candidates = [
        uvVenvPython,
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
