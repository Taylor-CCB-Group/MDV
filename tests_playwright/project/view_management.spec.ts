import test, { expect } from "@playwright/test";
import {
    createTemporaryProjectViaSyntheticAnndata,
    waitForProjectReady,
} from "../utils/projectFixtures";
import {
    createViewViaUi, getAllViews, getCurrentView, selectViewViaUi,
} from "../utils/helpers";

test.describe("View Management", () => {
    test.setTimeout(180_000);

    test("loads the default view in a synthetic AnnData project", async ({ page }) => {
        const projectHandle = await createTemporaryProjectViaSyntheticAnndata(page, {
            synthetic: {
                profile: "minimal",
                nCells: 200,
                nGenes: 12,
                force: true,
            },
        });

        try {
            await expect.poll(async () => await getCurrentView(page)).toBe("default");
            await expect.poll(async () => await getAllViews(page)).toContain("default");
        } finally {
            await projectHandle.cleanup();
        }
    });

    test("creates a new view and switches between it and default", async ({ page }) => {
        const projectHandle = await createTemporaryProjectViaSyntheticAnndata(page, {
            synthetic: {
                profile: "minimal",
                nCells: 200,
                nGenes: 12,
                force: true,
            },
        });

        try {
            const newViewName = `playwright_view_${Date.now()}`;
            const viewsBefore = await getAllViews(page);

            await createViewViaUi(page, newViewName);

            const viewsAfterCreate = await getAllViews(page);
            expect(viewsAfterCreate.length).toBe(viewsBefore.length + 1);
            expect(viewsAfterCreate).toContain(newViewName);

            await selectViewViaUi(page, "default");
            await selectViewViaUi(page, newViewName);
        } finally {
            await projectHandle.cleanup();
        }
    });

    test("deletes a created view and keeps the default view after reload", async ({ page }) => {
        const projectHandle = await createTemporaryProjectViaSyntheticAnndata(page, {
            synthetic: {
                profile: "minimal",
                nCells: 200,
                nGenes: 12,
                force: true,
            },
        });

        try {
            const newViewName = `delete_view_${Date.now()}`;

            await createViewViaUi(page, newViewName);
            const viewsBeforeDelete = await getAllViews(page);
            expect(viewsBeforeDelete).toContain(newViewName);

            await page.getByRole("button", { name: /delete current view/i }).click();
            const dialog = page.getByRole("dialog");
            await expect(dialog).toBeVisible();
            await dialog.getByRole("button", { name: "Yes", exact: true }).click();
            await expect(dialog).not.toBeVisible();

            await expect.poll(async () => await getCurrentView(page)).toBe("default");
            await expect.poll(async () => await getAllViews(page)).not.toContain(newViewName);

            await page.reload();
            await waitForProjectReady(page);

            await expect.poll(async () => await getCurrentView(page)).toBe("default");
            await expect.poll(async () => await getAllViews(page)).not.toContain(newViewName);
        } finally {
            await projectHandle.cleanup();
        }
    });
});
