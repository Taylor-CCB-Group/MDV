/**
 * `generate_synthetic_anndata_project.py` → folder under `~/mdv` → `/rescan_projects` → registered project id.
 */

import { expect, type Browser, type Page } from "@playwright/test";
import { execFile } from "node:child_process";
import fs from "node:fs/promises";
import os from "node:os";
import path from "node:path";
import { promisify } from "node:util";
import { getProjectUrl } from "../projectUrl";
import {
    deleteProjectViaApi,
    ensurePlaywrightMdvRoot,
    listProjectsViaApi,
    shouldDeleteProjectOnCleanup,
    triggerRescanProjects,
    waitForProjectReady,
} from "./core";
import { buildPythonPath, getPythonCandidates, REPO_ROOT } from "./pythonEnv";

const execFileAsync = promisify(execFile);

/** Mirrors `PROFILE_CHOICES` in `generate_synthetic_anndata_project.py`. */
export type SyntheticAnndataProfile = "minimal" | "realistic" | "large" | "memory-efficient";

/** CLI flags forwarded to `python -m mdvtools.tests.generate_synthetic_anndata_project`. */
export type SyntheticAnndataGeneratorCliArgs = {
    profile?: SyntheticAnndataProfile;
    nCells?: number;
    nGenes?: number;
    seed?: number;
    addMissing?: boolean;
    sparseDensity?: number;
    chunkData?: boolean;
    computeXUmap?: boolean;
    /** Adds synth_layer_a / synth_layer_b AnnData.layers → multiple rows_as_columns subgroups. */
    extraExpressionLayers?: boolean;
    force?: boolean;
};

export type SyntheticAnndataTemporaryProjectOptions = {
    deleteOnCleanup?: boolean;
    /** Must be a direct child of `~/mdv` when using the default scanner. */
    projectPath?: string;
    /** Used with `mdvRoot` when `projectPath` is omitted. */
    nameSegment?: string;
    mdvRoot?: string;
    synthetic?: SyntheticAnndataGeneratorCliArgs;
};

/** @deprecated Use {@link SyntheticAnndataTemporaryProjectOptions}. */
export type AnndataRescanTemporaryProjectOptions = SyntheticAnndataTemporaryProjectOptions;

export type SyntheticAnndataTemporaryProjectHandle = {
    cleanup: () => Promise<void>;
    projectId: string;
    projectUrl: string;
    projectPath: string;
    sourceUsed: "synthetic-anndata-filesystem-rescan";
};

export type SharedSyntheticAnndataSuiteHandle = SyntheticAnndataTemporaryProjectHandle & {
    openProjectPage: (page: Page) => Promise<void>;
};

/** @deprecated Use {@link SyntheticAnndataTemporaryProjectHandle}. */
export type AnndataRescanTemporaryProjectHandle = SyntheticAnndataTemporaryProjectHandle;

async function runSyntheticAnndataProjectCli(
    projectPath: string,
    synthetic: SyntheticAnndataGeneratorCliArgs,
) {
    const poetryProjectDir = path.join(REPO_ROOT, "python");
    const profile = synthetic.profile ?? "minimal";
    const nCells = synthetic.nCells ?? 1000;
    const nGenes = synthetic.nGenes ?? 50;
    const seed = synthetic.seed ?? 42;
    const sparseDensity = synthetic.sparseDensity ?? 0.1;
    const force = synthetic.force ?? true;

    const args = [
        "-m",
        "mdvtools.tests.generate_synthetic_anndata_project",
        "--profile",
        profile,
        "--n-cells",
        String(nCells),
        "--n-genes",
        String(nGenes),
        "--seed",
        String(seed),
        "--sparse-density",
        String(sparseDensity),
        "--output",
        projectPath,
    ];
    if (synthetic.addMissing) {
        args.push("--add-missing");
    }
    if (synthetic.chunkData) {
        args.push("--chunk-data");
    }
    if (synthetic.computeXUmap) {
        args.push("--compute-x-umap");
    }
    if (synthetic.extraExpressionLayers) {
        args.push("--extra-expression-layers");
    }
    if (force) {
        args.push("--force");
    }

    const errors: string[] = [];
    const pythonCandidates = await getPythonCandidates();
    for (const pythonPath of pythonCandidates) {
        try {
            await execFileAsync(pythonPath, args, {
                cwd: poetryProjectDir,
                env: {
                    ...process.env,
                    PYTHONPATH: buildPythonPath(),
                },
            });
            return;
        } catch (error) {
            const message = error instanceof Error ? error.message : String(error);
            errors.push(`${pythonPath}: ${message}`);
        }
    }

    throw new Error(
        `Unable to run mdvtools.tests.generate_synthetic_anndata_project.\n${errors.join("\n")}`,
    );
}

export async function createTemporaryProjectViaSyntheticAnndata(
    page: Page,
    options: SyntheticAnndataTemporaryProjectOptions = {},
): Promise<SyntheticAnndataTemporaryProjectHandle> {
    const deleteOnCleanup = shouldDeleteProjectOnCleanup(options.deleteOnCleanup);
    const mdvRoot = options.mdvRoot ?? (await ensurePlaywrightMdvRoot());

    const nameSegment =
        options.nameSegment ??
        `synth-anndata--playwright-rescan--minimal--${Date.now()}`;
    const projectPath = options.projectPath ?? path.join(mdvRoot, nameSegment);

    let projectId: string | undefined;

    try {
        const projectsBefore = await listProjectsViaApi(page.request);
        const projectIdsBefore = new Set(projectsBefore.map((project) => String(project.id)));

        await runSyntheticAnndataProjectCli(projectPath, options.synthetic ?? {});
        await triggerRescanProjects(page.request);

        const projectsAfter = await listProjectsViaApi(page.request);
        const createdProject =
            projectsAfter.find((project) => project.path === projectPath) ??
            projectsAfter.find((project) => project.name === nameSegment) ??
            projectsAfter.find((project) => !projectIdsBefore.has(String(project.id)));
        expect(createdProject).toBeTruthy();
        if (!createdProject) {
            throw new Error(
                "Synthetic AnnData project was generated but was not registered by /rescan_projects.",
            );
        }
        projectId = String(createdProject.id);

        const projectUrl = getProjectUrl(projectId);
        const response = await page.goto(projectUrl, { waitUntil: "domcontentloaded" });
        if (response && !response.ok()) {
            throw new Error(
                `Project page did not load: GET ${projectUrl} → ${response.status()} ${response.statusText()}. ` +
                    `If you see 500, check that browser context has the same baseURL as playwright.config (so /projects and rescan hit the real server). See server logs for the project id ${projectId}.`,
            );
        }
        await waitForProjectReady(page);

        return {
            cleanup: async () => {
                let cleanupError: unknown;
                if (deleteOnCleanup) {
                    try {
                        await deleteProjectViaApi(page.request, projectId);
                    } catch (error) {
                        cleanupError = error;
                    }
                }
                await fs.rm(projectPath, { recursive: true, force: true });
                if (cleanupError) {
                    throw cleanupError;
                }
            },
            projectId,
            projectUrl,
            projectPath,
            sourceUsed: "synthetic-anndata-filesystem-rescan",
        };
    } catch (error) {
        let cleanupError: unknown;
        if (deleteOnCleanup) {
            try {
                await deleteProjectViaApi(page.request, projectId);
            } catch (deleteError) {
                cleanupError = deleteError;
            }
        }
        await fs.rm(projectPath, { recursive: true, force: true }).catch(() => {});
        if (cleanupError) {
            throw cleanupError;
        }
        throw error;
    }
}

export async function createSharedSyntheticAnndataSuite(
    browser: Browser,
    options: SyntheticAnndataTemporaryProjectOptions = {},
): Promise<SharedSyntheticAnndataSuiteHandle> {
    const bootstrapPage = await browser.newPage();
    let projectHandle: SyntheticAnndataTemporaryProjectHandle | undefined;

    try {
        projectHandle = await createTemporaryProjectViaSyntheticAnndata(bootstrapPage, options);
        const handle = projectHandle;

        return {
            ...handle,
            cleanup: async () => {
                let cleanupError: unknown;
                try {
                    await handle.cleanup();
                } catch (error) {
                    cleanupError = error;
                }
                await bootstrapPage.close();
                if (cleanupError) {
                    throw cleanupError;
                }
            },
            openProjectPage: async (page: Page) => {
                await page.goto(handle.projectUrl, { waitUntil: "domcontentloaded" });
                await waitForProjectReady(page);
            },
        };
    } catch (error) {
        await bootstrapPage.close().catch(() => {});
        throw error;
    }
}

/** @deprecated Use {@link createTemporaryProjectViaSyntheticAnndata}. */
export const createTemporaryProjectViaAnndataRescan = createTemporaryProjectViaSyntheticAnndata;
