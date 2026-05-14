import { test, expect } from "@playwright/test";
import {
    gotoPath,
    mockApiRoot,
    mockExtensionConfig,
    mockProjects,
} from "../utils/routes";

test("import dialog opens when import permissions are enabled", async ({ page }) => {
    test.skip(
        (process.env.TEST_BASE_URL || "").includes("5055"),
        "Import-project card is hidden in production-style builds; run this spec against local Vite where import.meta.env.PROD is false.",
    );

    await mockApiRoot(page);
    await mockExtensionConfig(page);
    await mockProjects(page, []);

    await gotoPath(page);
    await page.waitForLoadState("networkidle");

    await expect(page.getByText("Import an existing project")).toBeVisible();
    await page.getByText("Import an existing project").click();

    const dialog = page.getByRole("dialog", { name: "Import Project" });
    await expect(dialog).toBeVisible();
    await expect(dialog.getByText(/drag and drop a file here/i)).toBeVisible();
    await expect(dialog.getByRole("button", { name: /close/i })).toBeVisible();

    await dialog.getByRole("button", { name: /close/i }).click();
    await expect(dialog).not.toBeVisible();
});
