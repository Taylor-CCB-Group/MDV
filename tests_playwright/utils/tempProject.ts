import { expect, type APIRequestContext, type Page } from "@playwright/test";
import { execFile } from "node:child_process";
import fs from "node:fs/promises";
import os from "node:os";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { promisify } from "node:util";
import { getProjectUrl } from "./projectUrl";

const execFileAsync = promisify(execFile);
const THIS_DIR = path.dirname(fileURLToPath(import.meta.url));
const REPO_ROOT = path.resolve(THIS_DIR, "..", "..");
const DEFAULT_PROJECT_NAME = "Playwright Temp Project";
const DEFAULT_FALLBACK_CSV = `sample_id,cell_type,score,cluster
s1,T-cell,0.91,c1
s2,B-cell,0.77,c2
s3,Myeloid,0.63,c1
s4,T-cell,0.88,c3
s5,NK,0.55,c2
`;

type MockProjectConfig = {
    nCells?: number;
    nGenes?: number;
};

export type TemporaryProjectOptions = {
    allowCsvFallback?: boolean;
    fallbackCsv?: string;
    mockConfig?: MockProjectConfig;
    projectName?: string;
};

export type TemporaryProjectHandle = {
    cleanup: () => Promise<void>;
    projectId: string;
    projectUrl: string;
    sourceUsed: "generated-mock" | "csv-fallback";
};

function getPythonCandidates() {
    const repoVenvPython = path.join(REPO_ROOT, "venv", "bin", "python");
    const candidates = [
        repoVenvPython,
        process.env.PYTHON,
        "python3",
        "python",
    ];
    return candidates.filter((candidate): candidate is string => Boolean(candidate));
}

function buildPythonPath() {
    const repoPythonPath = path.join(REPO_ROOT, "python");
    const existing = process.env.PYTHONPATH;
    return existing ? `${repoPythonPath}${path.delimiter}${existing}` : repoPythonPath;
}

async function generateMockProjectArchive(mockConfig: MockProjectConfig = {}) {
    const nCells = mockConfig.nCells ?? 8;
    const nGenes = mockConfig.nGenes ?? 12;
    const archivePath = path.join(
        os.tmpdir(),
        `mdv-playwright-${Date.now()}-${Math.random().toString(36).slice(2)}.mdv.zip`,
    );
    const pythonCode = [
        "from mdvtools.tests.test_project_factory import create_test_project_zip",
        `path = create_test_project_zip(r'''${archivePath}''', source='mock', n_cells=${nCells}, n_genes=${nGenes})`,
        "print(path)",
    ].join("\n");

    const errors: string[] = [];
    for (const candidate of getPythonCandidates()) {
        try {
            const { stdout } = await execFileAsync(candidate, ["-c", pythonCode], {
                cwd: REPO_ROOT,
                env: {
                    ...process.env,
                    PYTHONPATH: buildPythonPath(),
                },
            });
            if (stdout.trim()) {
                return archivePath;
            }
        } catch (error) {
            const message = error instanceof Error ? error.message : String(error);
            errors.push(`${candidate}: ${message}`);
        }
    }

    throw new Error(
        `Unable to generate a mock MDV project archive with the available Python environments.\n${errors.join("\n")}`,
    );
}

async function importProjectArchive(
    request: APIRequestContext,
    archivePath: string,
    projectName: string,
) {
    const fixtureBuffer = await fs.readFile(archivePath);
    const response = await request.post("/import_project", {
        multipart: {
            name: projectName,
            file: {
                name: path.basename(archivePath),
                mimeType: "application/zip",
                buffer: fixtureBuffer,
            },
        },
    });
    expect(response.ok()).toBe(true);
    const payload = await response.json();
    expect(payload?.status).toBe("success");
    expect(payload?.id).toBeTruthy();
    return String(payload.id);
}

async function createProjectViaApi(request: APIRequestContext) {
    const response = await request.post("/create_project");
    expect(response.ok()).toBe(true);
    const payload = await response.json();
    expect(payload?.id).toBeTruthy();
    return String(payload.id);
}

async function createFallbackCsvFile(csvContent: string) {
    const csvPath = path.join(
        os.tmpdir(),
        `mdv-playwright-${Date.now()}-${Math.random().toString(36).slice(2)}.csv`,
    );
    await fs.writeFile(csvPath, csvContent, "utf8");
    return csvPath;
}

async function uploadCsvToProject(page: Page, projectId: string, csvPath: string) {
    await page.goto(getProjectUrl(projectId));
    await page.waitForLoadState("networkidle");

    const fileInput = page.locator('input[type="file"]');
    await fileInput.waitFor({ state: "attached", timeout: 60_000 });
    await fileInput.setInputFiles(csvPath);

    await page.getByRole("button", { name: "Upload", exact: true }).waitFor({
        state: "visible",
        timeout: 30_000,
    });
    await page.getByRole("button", { name: "Upload", exact: true }).click();
    await page
        .getByRole("button", { name: /refresh page/i })
        .waitFor({ state: "visible", timeout: 120_000 });
    await page.getByRole("button", { name: /refresh page/i }).click();
}

export async function waitForProjectReady(page: Page) {
    await page.waitForLoadState("networkidle");
    await page.waitForFunction(() => Boolean((window as any).mdv?.chartManager?.viewManager));
    await expect(page.locator(".ciview-contentDiv").first()).toBeVisible({ timeout: 60_000 });
}

export async function deleteProjectViaApi(request: APIRequestContext, projectId: string | undefined) {
    if (!projectId) {
        return;
    }
    await request.delete(`/delete_project/${projectId}`);
}

export async function createTemporaryProject(
    page: Page,
    options: TemporaryProjectOptions = {},
): Promise<TemporaryProjectHandle> {
    const tempPaths: string[] = [];
    const projectName = options.projectName ?? DEFAULT_PROJECT_NAME;
    let projectId: string | undefined;
    let sourceUsed: "generated-mock" | "csv-fallback" = "generated-mock";

    try {
        try {
            const archivePath = await generateMockProjectArchive(options.mockConfig);
            tempPaths.push(archivePath);
            projectId = await importProjectArchive(page.request, archivePath, projectName);
        } catch (error) {
            if (options.allowCsvFallback === false) {
                throw error;
            }
            sourceUsed = "csv-fallback";
            const csvPath = await createFallbackCsvFile(options.fallbackCsv ?? DEFAULT_FALLBACK_CSV);
            tempPaths.push(csvPath);
            projectId = await createProjectViaApi(page.request);
            await uploadCsvToProject(page, projectId, csvPath);
        }

        const projectUrl = getProjectUrl(projectId);
        await page.goto(projectUrl);
        await waitForProjectReady(page);

        return {
            cleanup: async () => {
                await deleteProjectViaApi(page.request, projectId);
                await Promise.all(
                    tempPaths.map(async (tempPath) => {
                        await fs.rm(tempPath, { force: true });
                    }),
                );
            },
            projectId,
            projectUrl,
            sourceUsed,
        };
    } catch (error) {
        await deleteProjectViaApi(page.request, projectId);
        await Promise.all(
            tempPaths.map(async (tempPath) => {
                await fs.rm(tempPath, { force: true });
            }),
        );
        throw error;
    }
}
