/**
 * Run analysis: jobs submitted from the dialog run on the server and their output columns reach
 * the open page without a reload (ADR-0012, ADR-0009).
 *
 * The minimal synthetic project stores its expression matrix dense, so UMAP here also covers the
 * dense matrix tray.
 *
 * Run:
 *   TEST_BASE_URL=http://127.0.0.1:5170 pnpm run playwright-test-project tests_playwright/project/jobs/run_analysis.spec.ts --project=chromium --reporter=list
 */

import test, { expect, type Page } from "@playwright/test";
import {
    type SyntheticAnndataTemporaryProjectHandle,
    createTemporaryProjectViaSyntheticAnndata,
} from "../../utils/projectFixtures";

const SYNTHETIC = { profile: "minimal", nCells: 200, nGenes: 12, force: true } as const;
const JOB_TIMEOUT = 120_000;

/** Submit UMAP from the Run analysis dialog with the given output base name, then close the dialog. */
async function runUmap(page: Page, base: string, minDist?: string) {
    await page.getByRole("button", { name: "Run Analysis" }).click();
    const dialog = page.getByRole("dialog", { name: "Run analysis" });
    await dialog.getByRole("combobox", { name: "Tool" }).selectOption("umap");
    await dialog.getByRole("textbox", { name: "New column base name" }).fill(base);
    if (minDist !== undefined) {
        await dialog.getByRole("textbox", { name: "Minimum distance" }).fill(minDist);
    }
    await dialog.getByRole("button", { name: "Run", exact: true }).click();
    await expect(dialog.getByRole("row", { name: /UMAP/ }).first()).toBeVisible();
    await page.keyboard.press("Escape");
    await expect(dialog).toBeHidden();
}

/** The page's cells columns whose field starts with `base`, with the job that wrote each. */
async function jobColumns(page: Page, base: string) {
    return await page.evaluate((prefix): { field: string; job: string | undefined }[] => {
        const store = (window as any).mdv.chartManager.dsIndex.cells.dataStore;
        return store.columns
            .filter((c: any) => c.field.startsWith(prefix))
            .map((c: any) => ({ field: String(c.field), job: c.provenance?.job_id as string | undefined }));
    }, base);
}

test.describe("Run analysis", () => {
    test.setTimeout(180_000);

    test("running UMAP adds its columns to the open page", async ({ page }) => {
        let handle: SyntheticAnndataTemporaryProjectHandle | undefined;
        try {
            handle = await createTemporaryProjectViaSyntheticAnndata(page, {
                nameSegment: `jobs-umap--${Date.now()}`,
                synthetic: SYNTHETIC,
            });
            const base = `UMAP${Date.now()}`;

            await runUmap(page, base);

            // the dialog is closed, so this also checks the columns arrive without it
            await expect
                .poll(async () => (await jobColumns(page, base)).map((c) => c.field), { timeout: JOB_TIMEOUT })
                .toEqual([`${base}_1`, `${base}_2`]);
            const [first, second] = await jobColumns(page, base);
            expect(first.job).toBeTruthy();
            expect(second.job).toBe(first.job);

            await page.getByRole("button", { name: "Run Analysis" }).click();
            const row = page.getByRole("dialog", { name: "Run analysis" }).getByRole("row", { name: /UMAP/ });
            await expect(row.getByText("done")).toBeVisible();
            await expect(row.getByRole("cell").last()).toHaveText("");
        } finally {
            if (handle !== undefined) {
                await handle.cleanup();
            }
        }
    });

    test("re-running UMAP replaces its columns in place", async ({ page }) => {
        let handle: SyntheticAnndataTemporaryProjectHandle | undefined;
        try {
            handle = await createTemporaryProjectViaSyntheticAnndata(page, {
                nameSegment: `jobs-rerun--${Date.now()}`,
                synthetic: SYNTHETIC,
            });
            const base = `UMAP${Date.now()}`;
            const field = `${base}_1`;

            await runUmap(page, base);
            await expect.poll(async () => (await jobColumns(page, base)).length, { timeout: JOB_TIMEOUT }).toBe(2);
            const [{ job: firstJob }] = await jobColumns(page, base);

            // load the column's data and keep a reference to the column object, as a chart would
            const before = await page.evaluate(async (f) => {
                const mdv = (window as any).mdv;
                await new Promise<void>((resolve) => mdv.chartManager.loadColumnSet([f], "cells", resolve));
                const column = mdv.chartManager.dsIndex.cells.dataStore.columnIndex[f];
                (window as any).__rerunColumn = column;
                return Array.from(column.data.slice(0, 20) as Float32Array);
            }, field);

            await runUmap(page, base, "0.9");

            await expect
                .poll(async () => (await jobColumns(page, base)).map((c) => c.job), { timeout: JOB_TIMEOUT })
                .not.toContain(firstJob);
            const after = await page.evaluate((f) => {
                const column = (window as any).mdv.chartManager.dsIndex.cells.dataStore.columnIndex[f];
                return {
                    sameObject: column === (window as any).__rerunColumn,
                    values: Array.from(column.data.slice(0, 20) as Float32Array),
                };
            }, field);
            expect(after.sameObject).toBe(true);
            expect(after.values).not.toEqual(before);
            expect((await jobColumns(page, base)).map((c) => c.field)).toEqual([`${base}_1`, `${base}_2`]);
        } finally {
            if (handle !== undefined) {
                await handle.cleanup();
            }
        }
    });

    // ADR-0012: one JobManager per project, so a job belongs only to the project it was submitted to
    test("jobs in two open projects stay in their own project", async ({ page }) => {
        test.setTimeout(240_000);
        const pageB = await page.context().newPage();
        let handleA: SyntheticAnndataTemporaryProjectHandle | undefined;
        let handleB: SyntheticAnndataTemporaryProjectHandle | undefined;
        try {
            handleA = await createTemporaryProjectViaSyntheticAnndata(page, {
                nameSegment: `jobs-scope-a--${Date.now()}`,
                synthetic: SYNTHETIC,
            });
            handleB = await createTemporaryProjectViaSyntheticAnndata(pageB, {
                nameSegment: `jobs-scope-b--${Date.now()}`,
                synthetic: SYNTHETIC,
            });
            const a = { page, id: handleA.projectId, base: `UMAPa${Date.now()}` };
            const b = { page: pageB, id: handleB.projectId, base: `UMAPb${Date.now()}` };

            // both submitted before either finishes, so the one driver loop runs them together
            await runUmap(a.page, a.base);
            await runUmap(b.page, b.base);

            const jobIds: string[] = [];
            for (const p of [a, b]) {
                await expect
                    .poll(async () => (await jobColumns(p.page, p.base)).length, { timeout: JOB_TIMEOUT })
                    .toBe(2);
                const listed = await (await page.request.get(`/project/${p.id}/jobs`)).json();
                expect(listed).toHaveLength(1);
                expect((await jobColumns(p.page, p.base))[0].job).toBe(listed[0].job_id);
                jobIds.push(listed[0].job_id);
            }
            const [jobA, jobB] = jobIds;
            expect(jobA).not.toBe(jobB);

            // neither page gains the other project's columns
            expect(await jobColumns(a.page, b.base)).toEqual([]);
            expect(await jobColumns(b.page, a.base)).toEqual([]);

            // neither project's routes find the other project's job
            expect((await page.request.get(`/project/${a.id}/jobs/${jobB}`)).status()).toBe(404);
            expect((await page.request.get(`/project/${b.id}/jobs/${jobA}`)).status()).toBe(404);

            // each project's saved metadata holds only its own job's columns
            for (const [p, other, job] of [
                [a, b, jobA],
                [b, a, jobB],
            ] as const) {
                const configs = await (await page.request.get(`/project/${p.id}/get_configs`)).json();
                const cells = configs.datasources.find((d: { name: string }) => d.name === "cells");
                const fields = cells.columns.map((c: { field: string }) => c.field);
                expect(fields.filter((f: string) => f.startsWith(other.base))).toEqual([]);
                const own = cells.columns.filter((c: { field: string }) => c.field.startsWith(p.base));
                expect(own.map((c: { provenance?: { job_id?: string } }) => c.provenance?.job_id)).toEqual([job, job]);
            }
        } finally {
            if (handleB !== undefined) {
                await handleB.cleanup();
            }
            if (handleA !== undefined) {
                await handleA.cleanup();
            }
            await pageB.close();
        }
    });
});
