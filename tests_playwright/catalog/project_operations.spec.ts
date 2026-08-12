import { test, expect } from "@playwright/test";
import { gotoPath, mockApiRoot, mockExtensionConfig, mockProjects } from "../utils/routes";
import { project } from "../utils/data";

test.describe("Project operations", () => {
    test.beforeEach(async ({ page }) => {
        await mockApiRoot(page);
        await mockExtensionConfig(page);
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
        await expect(page.getByText("move this project to the recycle bin")).toBeVisible();
        await page.getByRole("button", { name: /move to recycle bin/i }).click();
        await expect(page.getByText("Alpha")).toHaveCount(0);
    });

    test("select multiple projects and move them to the recycle bin", async ({ page }) => {
        await mockProjects(page, [
            project({ id: 1, name: "Alpha" }),
            project({ id: 2, name: "Beta" }),
        ]);
        await page.reload();
        await page.route("**/projects/soft-delete", (route) =>
            route.fulfill({ json: { deletedProjectIds: [1, 2] } }),
        );

        await page.getByRole("button", { name: "Select" }).click();
        await page.getByLabel("select project Alpha").click();
        await page.getByLabel("select project Beta").click();
        await expect(page.getByText("2 selected")).toBeVisible();
        await page.getByRole("button", { name: "Delete selected" }).click();
        await page.getByRole("button", { name: "Move to recycle bin" }).click();

        await expect(page.getByText("Alpha")).toHaveCount(0);
        await expect(page.getByText("Beta")).toHaveCount(0);
    });

    test("preserves selection when switching between grid and list views", async ({ page }) => {
        await page.getByRole("button", { name: "Select" }).click();
        await page.getByLabel("select project Alpha").click();
        await page.getByRole("button", { name: /list view/i }).click();
        await expect(page.getByLabel("select project Alpha")).toBeChecked();
    });

    test("opens the recycle bin and permanently deletes all projects", async ({ page }) => {
        await page.route("**/projects/recycle-bin", async (route) => {
            if (route.request().method() === "DELETE") {
                await route.fulfill({ json: { deletedProjectIds: [2], failures: [] } });
                return;
            }
            await route.fulfill({
                json: [{ id: 2, name: "Beta", deletedTimestamp: "2025-01-01T00:00:00" }],
            });
        });

        await page.getByRole("button", { name: "Recycle bin" }).click();
        await expect(page.getByRole("dialog", { name: "Recycle bin" }).getByText("Beta")).toBeVisible();
        await page.getByRole("button", { name: "Delete all permanently" }).click();
        await page.getByRole("button", { name: "Permanently delete all" }).click();
        await expect(page.getByText("Your recycle bin is empty")).toBeVisible();
    });

    test("restores a project from the recycle bin", async ({ page }) => {
        await page.route("**/projects/recycle-bin/2/restore", (route) =>
            route.fulfill({ json: { restoredProjectId: 2 } }),
        );
        await page.route("**/projects/recycle-bin", (route) =>
            route.fulfill({
                json: [{ id: 2, name: "Beta", deletedTimestamp: "2025-01-01T00:00:00" }],
            }),
        );

        await page.getByRole("button", { name: "Recycle bin" }).click();
        await page.getByRole("button", { name: "Restore" }).click();
        await expect(page.getByRole("dialog", { name: "Restore project?" })).toBeVisible();
        await page.getByRole("button", { name: "Restore project" }).click();
        await expect(page.getByText("Your recycle bin is empty")).toBeVisible();
    });

    test("shows partial recycle bin purge failures", async ({ page }) => {
        await page.route("**/projects/recycle-bin", async (route) => {
            if (route.request().method() === "DELETE") {
                await route.fulfill({
                    json: { deletedProjectIds: [], failures: [{ projectId: 2, message: "busy" }] },
                });
                return;
            }
            await route.fulfill({ json: [{ id: 2, name: "Beta" }] });
        });

        await page.getByRole("button", { name: "Recycle bin" }).click();
        await page.getByRole("button", { name: "Delete all permanently" }).click();
        await page.getByRole("button", { name: "Permanently delete all" }).click();
        await expect(page.getByText(/2: busy/)).toBeVisible();
    });
});
