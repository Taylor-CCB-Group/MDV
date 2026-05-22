import test, { expect } from "@playwright/test";
import fs from "node:fs/promises";
import os from "node:os";
import path from "node:path";
import {
    createTemporaryProjectViaSyntheticAnndata,
    listProjectsViaApi,
} from "../utils/projectFixtures";

const LEAK_TEST_PREFIX = "intentional-cleanup-failure--synth-anndata";
const DEFAULT_SYNTH_PREFIX = "synth-anndata--playwright-rescan--minimal--";

// Test intentional failure for project cleanup
// Skipped by default
test.describe.skip("Manual project cleanup failure repro", () => {
    test.setTimeout(180_000);

    test("intentionally fails after project creation but should not leak the project", async ({
        page,
    }) => {
        const nameSegment = `${LEAK_TEST_PREFIX}--${Date.now()}`;
        const projectHandle = await createTemporaryProjectViaSyntheticAnndata(page, {
            nameSegment,
            synthetic: {
                profile: "minimal",
                nCells: 200,
                nGenes: 12,
                force: true,
            },
        });

        let cleanupCompleted = false;
        try {
            await expect(page.locator(".ciview-contentDiv").first()).toBeVisible();
            throw new Error("Intentional failure after synthetic project creation.");
        } finally {
            await projectHandle.cleanup();
            cleanupCompleted = true;

            const projects = await listProjectsViaApi(page.request);
            const stillVisible = projects.some(
                (project) =>
                    project.name === nameSegment ||
                    String(project.id) === projectHandle.projectId,
            );
            expect(stillVisible).toBe(false);
            expect(
                projects.filter((project) => project.name?.startsWith(DEFAULT_SYNTH_PREFIX)),
            ).toEqual([]);

            const generatedPath = path.join(os.homedir(), "mdv", nameSegment);
            await expect(fs.stat(generatedPath).then(
                () => true,
                () => false,
            )).resolves.toBe(false);
            expect(await findGeneratedFolders(LEAK_TEST_PREFIX)).toEqual([]);
            expect(await findGeneratedFolders(DEFAULT_SYNTH_PREFIX)).toEqual([]);
            expect(cleanupCompleted).toBe(true);
        }
    });
});

async function findGeneratedFolders(prefix: string) {
    const mdvRoot = path.join(os.homedir(), "mdv");
    const entries = await fs.readdir(mdvRoot, { withFileTypes: true });
    return entries
        .filter((entry) => entry.isDirectory() && entry.name.startsWith(prefix))
        .map((entry) => entry.name);
}
