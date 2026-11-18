import fs from "node:fs/promises";
import os from "node:os";
import path from "node:path";
import { type APIRequestContext, expect, type Page } from "@playwright/test";

const DELETE_PROJECT_RETRIES = 3;
const DELETE_PROJECT_BACKOFF_MS = 1_000;

export function shouldKeepProjectsOnCleanup() {
    return process.env.PLAYWRIGHT_KEEP_PROJECTS === "1";
}

export function shouldDeleteProjectOnCleanup(explicit?: boolean) {
    if (explicit !== undefined) {
        return explicit;
    }
    return !shouldKeepProjectsOnCleanup();
}

export async function ensurePlaywrightMdvRoot() {
    const mdvRoot = path.join(os.homedir(), "mdv");
    await fs.mkdir(mdvRoot, { recursive: true });
    return mdvRoot;
}

/**
 * Generate a project somewhere the backend can't see it, then move it to `projectPath` in one step.
 *
 * Docker Desktop's file sharing takes a file's size when the host creates it and doesn't notice
 * it growing afterwards, so a project written straight into the shared `~/mdv` (h5 files are
 * created small and grown in place) reads as truncated in the container and `/rescan_projects`
 * fails to register it. The staging folder sits beside the shared folder so the rename stays on
 * one filesystem.
 */
export async function generateProjectOffShare(
    projectPath: string,
    generate: (stagingProjectPath: string) => Promise<void>,
) {
    const stagingRoot = await fs.mkdtemp(
        path.join(path.dirname(path.dirname(projectPath)), ".mdv-playwright-staging-"),
    );
    try {
        const stagingProjectPath = path.join(stagingRoot, path.basename(projectPath));
        await generate(stagingProjectPath);
        await fs.rm(projectPath, { recursive: true, force: true });
        await fs.rename(stagingProjectPath, projectPath);
    } finally {
        await fs.rm(stagingRoot, { recursive: true, force: true });
    }
}

async function sleep(ms: number) {
    await new Promise((resolve) => setTimeout(resolve, ms));
}

export async function waitForProjectReady(page: Page) {
    // Vite dev pages can keep HMR/network channels open, so networkidle may never resolve.
    await page.waitForLoadState("domcontentloaded");
    await page.waitForFunction(() => Boolean((window as any).mdv?.chartManager?.viewManager));
    await expect(page.locator(".ciview-contentDiv").first()).toBeVisible({ timeout: 60_000 });
    const forbiddenDialog = page.getByRole("dialog").filter({
        hasText: "You do not have permission to save this project.",
    });
    if (await forbiddenDialog.isVisible({ timeout: 1_000 }).catch(() => false)) {
        await forbiddenDialog.getByRole("button", { name: "close", exact: true }).click();
        await expect(forbiddenDialog).not.toBeVisible();
    }
}

export async function deleteProjectViaApi(request: APIRequestContext, projectId: string | undefined) {
    if (!projectId) {
        return;
    }
    let lastError: unknown;
    for (let attempt = 1; attempt <= DELETE_PROJECT_RETRIES; attempt += 1) {
        try {
            const response = await request.delete(`/delete_project/${projectId}`, {
                timeout: 15_000,
            });
            if (!response.ok()) {
                throw new Error(`DELETE /delete_project/${projectId} returned ${response.status()}`);
            }

            const projects = await listProjectsViaApi(request);
            const stillExists = projects.some((project) => String(project.id) === projectId);
            if (!stillExists) {
                return;
            }

            throw new Error(`project ${projectId} still present in /projects after delete`);
        } catch (error) {
            lastError = error;
            if (attempt < DELETE_PROJECT_RETRIES) {
                await sleep(DELETE_PROJECT_BACKOFF_MS * attempt);
            }
        }
    }

    const message = lastError instanceof Error ? lastError.message : String(lastError);
    throw new Error(`Failed to clean up Playwright project ${projectId}: ${message}`);
}

export type ProjectListSummary = {
    id: number | string;
    name?: string;
    path?: string;
};

export async function listProjectsViaApi(request: APIRequestContext): Promise<ProjectListSummary[]> {
    const response = await request.get("/projects");
    expect(response.ok()).toBe(true);
    const payload = await response.json();
    return Array.isArray(payload) ? (payload as ProjectListSummary[]) : [];
}

export async function triggerRescanProjects(request: APIRequestContext) {
    await request.get("/rescan_projects", {
        maxRedirects: 0,
        failOnStatusCode: false,
    });
}

function isProjectInFolder(project: ProjectListSummary, projectPath: string) {
    const folder = path.basename(projectPath);
    return project.name === folder || path.basename(project.path ?? "") === folder;
}

/**
 * Fail before generating into a folder the backend already has a project for: replacing its files
 * leaves that row pointing at them, so the rescan never registers a new project and cleanup, which
 * deletes by the new id, would leave the row behind.
 */
export function assertProjectPathUnregistered(projects: ProjectListSummary[], projectPath: string) {
    const existing = projects.find((project) => isProjectInFolder(project, projectPath));
    if (existing) {
        throw new Error(
            `Project ${existing.id} is already registered for ${projectPath}; use a fresh nameSegment / projectPath.`,
        );
    }
}

/**
 * Rescan until the project in `projectPath` is registered under a new id. The backend's view of a
 * Docker Desktop shared folder can lag the host, so a single rescan may not see a folder that was
 * just moved in (and may still see folders that were just deleted - hence matching on the folder
 * name rather than taking whichever project is new).
 */
export async function rescanUntilRegistered(
    request: APIRequestContext,
    projectPath: string,
    projectIdsBefore: Set<string>,
    timeoutMs = 30_000,
): Promise<ProjectListSummary | undefined> {
    const deadline = Date.now() + timeoutMs;
    for (;;) {
        await triggerRescanProjects(request);
        const created = (await listProjectsViaApi(request)).find(
            (project) => !projectIdsBefore.has(String(project.id)) && isProjectInFolder(project, projectPath),
        );
        if (created || Date.now() > deadline) {
            return created;
        }
        await sleep(1_000);
    }
}
