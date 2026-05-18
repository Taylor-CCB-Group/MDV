import test, { expect } from "@playwright/test";
import { createTemporaryProjectViaSyntheticAnndata } from "../../utils/projectFixtures";
import {
    addChartViaUi,
    expectChartPanelHasNoError,
    getNumericColumnNamesForChart,
    getSelectionToolbar,
    patchScatterChartConfig,
    waitForChartByTitle,
} from "../../utils/helpers";

test.describe("Density grid selection overlay", () => {
    test.setTimeout(180_000);

    test("keeps shape-drawing tools in grid mode and restores them in overlay mode", async ({ page }) => {
        const projectHandle = await createTemporaryProjectViaSyntheticAnndata(page, {
            synthetic: {
                profile: "minimal",
                nCells: 200,
                nGenes: 12,
                force: true,
            },
        });

        try {
            const title = `2D Scatter Density Grid ${Date.now()}`;
            await addChartViaUi(page, "2D Scatter Plot", title);
            await waitForChartByTitle(page, title);
            await expectChartPanelHasNoError(page, title);

            const densityFields = await getNumericColumnNamesForChart(page, title, 2);
            expect(densityFields.length).toBeGreaterThan(0);

            await patchScatterChartConfig(page, title, {
                type: "DeckContourScatter",
                densityFields,
                density_mode: "overlay",
            });

            const toolbar = getSelectionToolbar(page, title);
            await expect(toolbar.getByRole("button", { name: "Rectangle", exact: true })).toBeVisible();
            await expect(toolbar.getByRole("button", { name: "Pan", exact: true })).toBeVisible();
            await expect(toolbar.getByRole("button", { name: "Density grid view", exact: true })).toBeVisible();

            await toolbar.getByRole("button", { name: "Density grid view", exact: true }).click();

            await expect(toolbar.getByRole("button", { name: "Rectangle", exact: true })).toBeVisible();
            await expect(toolbar.getByRole("button", { name: "Polygon", exact: true })).toBeVisible();
            await expect(toolbar.getByRole("button", { name: "Freehand", exact: true })).toBeVisible();
            await expect(toolbar.getByRole("button", { name: "Modify", exact: true })).toBeVisible();
            await expect(toolbar.getByRole("button", { name: "Pan", exact: true })).toBeVisible();
            await expect(toolbar.getByRole("button", { name: "Single scatter view", exact: true })).toBeVisible();
            await expect(toolbar.getByRole("button", { name: "Manage Gates", exact: true })).toBeVisible();

            const panel = page.locator(".ciview-chart-panel").filter({ hasText: title }).first();
            await expect(panel.getByText("Choose density fields to build the grid.")).toHaveCount(0);
            await expect(panel.locator(".deckgl-overlay").first()).toBeVisible({ timeout: 30_000 });

            await toolbar.getByRole("button", { name: "Single scatter view", exact: true }).click();
            await expect(toolbar.getByRole("button", { name: "Rectangle", exact: true })).toBeVisible();
            await expect(toolbar.getByRole("button", { name: "Manage Gates", exact: true })).toBeVisible();
        } finally {
            await projectHandle.cleanup();
        }
    });
});
