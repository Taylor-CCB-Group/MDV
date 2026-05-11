import { test, expect } from "@playwright/test";
import {
    gotoPath,
    mockApiRoot,
    mockExtensionConfig,
    mockProjects,
} from "../utils/routes";

test("create project navigates to project view", async ({ page }) => {
    test.skip(
        (process.env.TEST_BASE_URL || "").includes("5055"),
        "Create-project card is hidden in production-style builds; run this spec against local Vite where import.meta.env.PROD is false.",
    );

    await mockApiRoot(page);
    await mockExtensionConfig(page);
    await mockProjects(page, []);

    await page.route("**/create_project", (route) =>
        route.fulfill({ json: { id: "123", name: "New Project" } }),
    );

    await gotoPath(page);

    await expect(page.getByText("Create new project")).toBeVisible();
    await page.getByText("Create new project").click();

    await page.waitForURL("**/project/123");
    expect(page.url()).toContain("/project/123");
});
