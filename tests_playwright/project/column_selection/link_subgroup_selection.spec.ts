/**
 * Rows-as-columns subgroup coverage for multilayer mock projects.
 *
 * Persistence and metadata checks run against the default backend (localhost:5055).
 * Subgroup tab UI checks need the worktree Vite frontend, e.g.:
 *   TEST_BASE_URL=http://127.0.0.1:5173 pnpm run playwright-test-project tests_playwright/project/column_selection/link_subgroup_selection.spec.ts --project=chromium
 */
import test, { expect } from "@playwright/test";
import {
    createTemporaryProjectViaSyntheticAnndata,
    waitForProjectReady,
} from "../../utils/projectFixtures";
import {
    addDotPlotChartViaApi,
    closeTopDialog,
    expectChartPanelHasNoError,
    getDotPlotXAxisFields,
    getFirstLinkedGeneLabel,
    getLinkColumnSettingRow,
    getRowsAsColumnsSubgroupKeys,
    openChartSettingsDialog,
    saveCurrentView,
    searchChartSettings,
    selectLinkedGeneInChartSettings,
    setDotPlotXAxisLinkField,
    isWorktreeFrontendForPlaywright,
    waitForChartByTitle,
    waitForViewUnsavedState,
} from "../../utils/helpers";

const MULTILAYER_SYNTHETIC = {
    profile: "minimal" as const,
    nCells: 80,
    nGenes: 40,
    seed: 0,
    extraExpressionLayers: true,
    force: true,
};

test.describe("Link subgroup column selection", () => {
    test.setTimeout(180_000);

    test("project exposes multiple rows_as_columns expression subgroups", async ({ page }) => {
        const projectHandle = await createTemporaryProjectViaSyntheticAnndata(page, {
            synthetic: MULTILAYER_SYNTHETIC,
        });

        try {
            const subgroups = await getRowsAsColumnsSubgroupKeys(page);
            const keys = subgroups.map((sg) => sg.key);
            expect(keys).toContain("gs");
            expect(keys).toContain("synth_layer_a");
            expect(keys).toContain("synth_layer_b");
            expect(keys.length).toBeGreaterThanOrEqual(3);
        } finally {
            await projectHandle.cleanup();
        }
    });

    test("shows subgroup layer tabs in the link column picker when the worktree UI is served", async ({
        page,
    }) => {
        const projectHandle = await createTemporaryProjectViaSyntheticAnndata(page, {
            synthetic: MULTILAYER_SYNTHETIC,
        });

        try {
            const subgroups = await getRowsAsColumnsSubgroupKeys(page);
            test.skip(
                !isWorktreeFrontendForPlaywright(),
                "Subgroup tabs are only in the worktree frontend; set TEST_BASE_URL to the Vite dev server.",
            );

            const title = `Dot Plot Subgroup Tabs ${Date.now()}`;
            await addDotPlotChartViaApi(page, title);
            await waitForChartByTitle(page, title);
            await openChartSettingsDialog(page, title);
            await searchChartSettings(page, "Fields on x axis");
            const { row, subgroupTablist } = getLinkColumnSettingRow(page, "Fields on x axis");
            await row.getByRole("tab", { name: "link", exact: true }).click();

            for (const { label } of subgroups) {
                await expect(subgroupTablist.getByRole("tab", { name: label, exact: true })).toBeVisible();
            }
            await closeTopDialog(page);
        } finally {
            await projectHandle.cleanup();
        }
    });

    test("persists dot plots using different link subgroups after save and reload", async ({ page }) => {
        const projectHandle = await createTemporaryProjectViaSyntheticAnndata(page, {
            synthetic: MULTILAYER_SYNTHETIC,
        });

        try {
            const geneLabel = await getFirstLinkedGeneLabel(page);
            const geneIndex = 0;
            const titleGs = `Dot Plot GS ${Date.now()}`;
            const titleLayerA = `Dot Plot Layer A ${Date.now() + 1}`;

            await addDotPlotChartViaApi(page, titleGs);
            await waitForChartByTitle(page, titleGs);

            const useSubgroupTabs = isWorktreeFrontendForPlaywright();
            if (useSubgroupTabs) {
                await openChartSettingsDialog(page, titleGs);
                await selectLinkedGeneInChartSettings(page, {
                    paramLabel: "Fields on x axis",
                    subgroupTabLabel: "Gene Scores",
                    geneLabel,
                });
                await closeTopDialog(page);
            } else {
                await setDotPlotXAxisLinkField(page, titleGs, "gs", geneLabel, geneIndex);
            }
            await waitForViewUnsavedState(page, true);
            await expectChartPanelHasNoError(page, titleGs);

            await addDotPlotChartViaApi(page, titleLayerA);
            await waitForChartByTitle(page, titleLayerA);
            if (useSubgroupTabs) {
                await openChartSettingsDialog(page, titleLayerA);
                await selectLinkedGeneInChartSettings(page, {
                    paramLabel: "Fields on x axis",
                    subgroupTabLabel: "synth_layer_a",
                    geneLabel,
                });
                await closeTopDialog(page);
            } else {
                await setDotPlotXAxisLinkField(page, titleLayerA, "synth_layer_a", geneLabel, geneIndex);
            }
            await waitForViewUnsavedState(page, true);
            await expectChartPanelHasNoError(page, titleLayerA);

            const gsFieldsBeforeSave = await getDotPlotXAxisFields(page, titleGs);
            const layerAFieldsBeforeSave = await getDotPlotXAxisFields(page, titleLayerA);
            expect(gsFieldsBeforeSave[0]).toMatch(/^gs\|/);
            expect(layerAFieldsBeforeSave[0]).toMatch(/^synth_layer_a\|/);

            await saveCurrentView(page);

            await page.reload();
            await waitForProjectReady(page);

            await expectChartPanelHasNoError(page, titleGs);
            await expectChartPanelHasNoError(page, titleLayerA);

            const gsFieldsAfterReload = await getDotPlotXAxisFields(page, titleGs);
            const layerAFieldsAfterReload = await getDotPlotXAxisFields(page, titleLayerA);
            expect(gsFieldsAfterReload).toEqual(gsFieldsBeforeSave);
            expect(layerAFieldsAfterReload).toEqual(layerAFieldsBeforeSave);
            expect(gsFieldsAfterReload[0]).toMatch(/^gs\|/);
            expect(layerAFieldsAfterReload[0]).toMatch(/^synth_layer_a\|/);
        } finally {
            await projectHandle.cleanup();
        }
    });
});
