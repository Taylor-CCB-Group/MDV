/**
 * `generate_synthetic_spatial_project.py` → folder under `~/mdv` → `/rescan_projects` → registered project id.
 */

import { expect, type Page } from "@playwright/test";
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

/** CLI flags forwarded to `python -m mdvtools.tests.generate_synthetic_spatial_project`. */
export type SyntheticSpatialGeneratorCliArgs = {
    profile?: string;
    nCells?: string | number;
    nGenes?: string | number;
    force?: boolean;
};

export type SyntheticSpatialTemporaryProjectOptions = {
    deleteOnCleanup?: boolean;
    /** When set, writes the project here (must be a direct child of `mdvRoot` / `~/mdv`). */
    projectPath?: string;
    /** Used with `mdvRoot` when `projectPath` is omitted. */
    nameSegment?: string;
    mdvRoot?: string;
    synthetic?: SyntheticSpatialGeneratorCliArgs;
};

/** @deprecated Use {@link SyntheticSpatialTemporaryProjectOptions}. */
export type RescanTemporaryProjectOptions = SyntheticSpatialTemporaryProjectOptions;

export type SyntheticSpatialTemporaryProjectHandle = {
    cleanup: () => Promise<void>;
    projectId: string;
    projectUrl: string;
    projectPath: string;
    sourceUsed: "synthetic-filesystem-rescan";
};

/** @deprecated Use {@link SyntheticSpatialTemporaryProjectHandle}. */
export type RescanTemporaryProjectHandle = SyntheticSpatialTemporaryProjectHandle;

async function runSyntheticSpatialProjectCli(
    projectPath: string,
    synthetic: SyntheticSpatialGeneratorCliArgs,
) {
    const poetryProjectDir = path.join(REPO_ROOT, "python");
    const profile = synthetic.profile ?? "scatter-table";
    const nCells = synthetic.nCells ?? "200";
    const nGenes = synthetic.nGenes ?? "12";
    const force = synthetic.force ?? true;

    const args = [
        "-m",
        "mdvtools.tests.generate_synthetic_spatial_project",
        "--profile",
        profile,
        "--n-cells",
        String(nCells),
        "--n-genes",
        String(nGenes),
        "--output",
        projectPath,
    ];
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
        `Unable to run mdvtools.tests.generate_synthetic_spatial_project.\n${errors.join("\n")}`,
    );
}

export async function createTemporaryProjectViaSyntheticSpatial(
    page: Page,
    options: SyntheticSpatialTemporaryProjectOptions = {},
): Promise<SyntheticSpatialTemporaryProjectHandle> {
    const deleteOnCleanup = shouldDeleteProjectOnCleanup(options.deleteOnCleanup);
    const mdvRoot = options.mdvRoot ?? (await ensurePlaywrightMdvRoot());

    const nameSegment =
        options.nameSegment ??
        `synth-spatial--playwright-rescan--scatter-table--${Date.now()}`;
    const projectPath = options.projectPath ?? path.join(mdvRoot, nameSegment);

    let projectId: string | undefined;

    try {
        const projectsBefore = await listProjectsViaApi(page.request);
        const projectIdsBefore = new Set(projectsBefore.map((project) => String(project.id)));

        await runSyntheticSpatialProjectCli(projectPath, options.synthetic ?? {});
        await triggerRescanProjects(page.request);

        const projectsAfter = await listProjectsViaApi(page.request);
        const createdProject =
            projectsAfter.find((project) => project.name === nameSegment) ??
            projectsAfter.find((project) => project.path === projectPath) ??
            projectsAfter.find((project) => !projectIdsBefore.has(String(project.id)));
        expect(createdProject).toBeTruthy();
        if (!createdProject) {
            throw new Error(
                "Synthetic project was generated but was not registered by /rescan_projects.",
            );
        }
        projectId = String(createdProject.id);

        const projectUrl = getProjectUrl(projectId);
        await page.goto(projectUrl);
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
            sourceUsed: "synthetic-filesystem-rescan",
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

/** @deprecated Use {@link createTemporaryProjectViaSyntheticSpatial}. */
export const createTemporaryProjectViaRescan = createTemporaryProjectViaSyntheticSpatial;
