import { expect, type APIRequestContext, type Page } from "@playwright/test";
import os from "node:os";
import path from "node:path";
import fs from "node:fs/promises";

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
