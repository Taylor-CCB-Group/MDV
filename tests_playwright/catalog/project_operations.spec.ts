import { test, expect } from "@playwright/test";
import { mockApiRoot, mockProjects, gotoPath } from "../utils/routes";
import { project } from "../utils/data";

test.describe("Project operations", () => {
    test.beforeEach(async ({ page }) => {
        await mockApiRoot(page);
        await mockProjects(page, [
            project({ id: "p1", name: "Alpha" }),
        ]);
        await gotoPath(page);
    });

    test("info project", async ({ page }) => {
        await page.getByTestId("project_menu_p1").click();
        await page.getByTestId("project_info_p1").click();
        await expect(page.getByText("Project Information")).toBeVisible();
        await expect(page.getByText("Readme")).toBeVisible();
        await page.getByRole("button", { name: /close/i }).click();
        await expect(page.getByText("Project Information")).not.toBeVisible();
    });

    test("rename project", async ({ page }) => {
        await page.route("**/projects/p1/rename", (r) => r.fulfill({ json: { ok: true } }));
        await page.getByTestId("project_menu_p1").click();
        await page.getByTestId("project_rename_p1").click();
        await page.getByLabel("Project Name").fill("Alpha Renamed");
        await page.getByRole("button", { name: /save changes/i }).click();
        await expect(page.getByText("Alpha Renamed")).toBeVisible();
    });

    test("delete project", async ({ page }) => {
        await page.route("**/delete_project/p1", (r) => r.fulfill({ json: { ok: true } }));
        await page.getByTestId("project_menu_p1").click();
        await page.getByTestId("project_delete_p1").click();
        await page.getByRole("button", { name: /delete project/i }).click();
        await expect(page.getByText("Alpha")).toHaveCount(0);
    });

    test("export project", async ({ page }) => {
        
      const downloadPromise = page.waitForEvent("download");

      await page.route("**/export_project/p1", (r) =>
            r.fulfill({ body: "ZIP", headers: { "Content-Type": "application/zip" } }),
        );
        await page.getByTestId("project_menu_p1").click();
        await page.getByTestId("project_export_p1").click();
        const dl = await downloadPromise;
        expect(dl.suggestedFilename()).toContain(".zip");
        expect(dl.suggestedFilename()).toEqual("Alpha.mdv.zip");
    });
});
