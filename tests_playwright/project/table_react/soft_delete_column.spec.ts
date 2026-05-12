import test, { expect, type Page } from "@playwright/test";
import {
    createTemporaryProjectViaSyntheticAnndata,
    type TemporaryProjectHandle,
    waitForProjectReady,
} from "../../utils/projectFixtures";

type ProjectSetupInfo = {
    cloneSourceName: string;
    dataSourceName: string;
    tableChartId: string;
};

async function setupReactTableProject(page: Page): Promise<ProjectSetupInfo> {
    return await page.evaluate(async () => {
        const chartManager = (window as any).mdv.chartManager;
        const dataSourceName = Object.keys(chartManager.dsIndex)[0];
        if (!dataSourceName) {
            throw new Error("No datasource available in project");
        }

        const dataStore = chartManager.dsIndex[dataSourceName].dataStore;
        const visibleColumns = dataStore.columns.filter(
            (column: any) => !column.subgroup && column.field !== "__index__",
        );
        if (visibleColumns.length < 3) {
            throw new Error("Need at least three columns to create the React table");
        }

        const cloneSource =
            visibleColumns.find((column: any) =>
                ["text", "text16", "unique", "multitext"].includes(column.datatype),
            ) ?? visibleColumns[0];
        const tableFields = visibleColumns.slice(0, 3).map((column: any) => column.field);
        const tableChartId = `soft-delete-table-${Date.now()}`;

        await chartManager.addChart(
            dataSourceName,
            {
                id: tableChartId,
                title: "Soft Delete Table",
                legend: "",
                type: "table_chart_react",
                param: tableFields,
                size: [900, 420],
                position: [20, 20],
            },
            true,
        );

        return {
            dataSourceName,
            tableChartId,
            cloneSourceName: cloneSource.name,
        };
    });
}

async function addClonedEditableColumn(
    page: Page,
    chartId: string,
    columnName: string,
    cloneSourceName: string,
) {
    await page.evaluate((targetChartId) => {
        const chart = Object.values((window as any).mdv.chartManager.charts)
            .map((entry: any) => entry.chart)
            .find((entry: any) => entry.config.id === targetChartId);
        if (!chart) {
            throw new Error(`React table chart ${targetChartId} not found`);
        }
        chart.openAddColumnDialog();
    }, chartId);

    const dialog = page.getByRole("dialog", { name: "Add Column" });
    await expect(dialog).toBeVisible();

    await dialog.getByLabel("Column Name").fill(columnName);
    await dialog.getByRole("checkbox", { name: "Clone existing column" }).check();

    const cloneInput = dialog.getByRole("combobox", { name: "Column to clone" });
    await cloneInput.click();
    await cloneInput.fill(cloneSourceName);
    await page
        .locator('[role="option"]')
        .filter({ hasText: cloneSourceName })
        .first()
        .click();

    await dialog.getByRole("button", { name: "Add", exact: true }).click();
    await expect(dialog).not.toBeVisible();

    await page.waitForFunction((targetColumnName) => {
        const chartManager = (window as any).mdv.chartManager;
        const dataSourceName = Object.keys(chartManager.dsIndex)[0];
        const dataStore = chartManager.dsIndex[dataSourceName].dataStore;
        const tableChart = Object.values(chartManager.charts)
            .map((entry: any) => entry.chart)
            .find((chart: any) => chart.config.title === "Soft Delete Table");
        return Boolean(
            dataStore.columnIndex[targetColumnName] &&
                Array.isArray(tableChart?.config?.param) &&
                tableChart.config.param.includes(targetColumnName),
        );
    }, columnName);
}

async function addDependentRowChart(
    page: Page,
    dataSourceName: string,
    columnName: string,
    title: string,
) {
    const rowChartId = await page.evaluate(
        async ({ dsName, field, chartTitle }) => {
            const chartManager = (window as any).mdv.chartManager;
            const rowChartId = `soft-delete-dependency-${Date.now()}`;
            await chartManager.addChart(
                dsName,
                {
                    id: rowChartId,
                    title: chartTitle,
                    legend: "",
                    type: "row_chart",
                    param: [field],
                    size: [420, 300],
                    position: [950, 20],
                },
                true,
            );
            return rowChartId;
        },
        { dsName: dataSourceName, field: columnName, chartTitle: title },
    );

    await page.waitForFunction((targetRowChartId) => {
        return Object.values((window as any).mdv.chartManager.charts).some(
            (entry: any) => entry.chart.config.id === targetRowChartId,
        );
    }, rowChartId);

    return rowChartId;
}
async function injectSavedViewImpact(
    page: Page,
    dataSourceName: string,
    columnName: string,
    savedViewName: string,
    chartTitle: string,
) {
    await page.evaluate(
        ({ targetDataSourceName, targetColumnName, targetViewName, targetChartTitle }) => {
            const chartManager = (window as any).mdv.chartManager;
            const viewManager = chartManager.viewManager;
            const allViews = Array.isArray(viewManager.all_views) ? viewManager.all_views : [];
            if (!allViews.includes(targetViewName)) {
                viewManager.all_views = [...allViews, targetViewName];
            }
            const originalViewLoader = chartManager.viewLoader?.bind(chartManager);
            chartManager.viewLoader = async (viewName: string) => {
                if (viewName === targetViewName) {
                    return {
                        dataSources: {
                            [targetDataSourceName]: {},
                        },
                        initialCharts: {
                            [targetDataSourceName]: [
                                {
                                    id: `saved-view-impact-${targetColumnName}`,
                                    title: targetChartTitle,
                                    legend: "",
                                    type: "row_chart",
                                    param: [targetColumnName],
                                    size: [420, 300],
                                    position: [950, 20],
                                },
                            ],
                        },
                    };
                }
                if (!originalViewLoader) {
                    return null;
                }
                return await originalViewLoader(viewName);
            };
        },
        {
            targetDataSourceName: dataSourceName,
            targetColumnName: columnName,
            targetViewName: savedViewName,
            targetChartTitle: chartTitle,
        },
    );
}

async function requestColumnDeletion(page: Page, chartId: string, columnName: string) {
    await page.evaluate(({ targetChartId, targetColumnName }) => {
        const chart = Object.values((window as any).mdv.chartManager.charts)
            .map((entry: any) => entry.chart)
            .find((entry: any) => entry.config.id === targetChartId);
        if (!chart?.gridRef?.current?.slickGrid) {
            throw new Error(`Grid not ready for chart ${targetChartId}`);
        }
        const slickGrid = chart.gridRef.current.slickGrid;
        const pubSub = slickGrid.getPubSubService?.();
        if (!pubSub?.publish) {
            throw new Error("SlickGrid PubSub service is unavailable");
        }
        const column =
            slickGrid.getColumns().find((entry: any) => entry.field === targetColumnName) ?? {
                field: targetColumnName,
                id: targetColumnName,
            };
        pubSub.publish("onHeaderMenuCommand", {
            column,
            command: "remove-column",
        });
    }, { targetChartId: chartId, targetColumnName: columnName });
}

async function saveCurrentView(page: Page) {
    await page.evaluate(async () => {
        const viewManager = (window as any).mdv.chartManager.viewManager;
        await viewManager.saveView();
    });
    await page.waitForFunction(() => {
        const viewManager = (window as any).mdv.chartManager.viewManager;
        return viewManager && viewManager.hasUnsavedChanges?.() === false;
    });
}

async function createSavedEditableColumn(
    page: Page,
    chartId: string,
    columnName: string,
    cloneSourceName: string,
) {
    await addClonedEditableColumn(page, chartId, columnName, cloneSourceName);
    await saveCurrentView(page);
    await page.reload();
    await waitForProjectReady(page);
    await page.waitForFunction((targetColumnName) => {
        const chartManager = (window as any).mdv.chartManager;
        const dataSourceName = Object.keys(chartManager.dsIndex)[0];
        const dataStore = chartManager.dsIndex[dataSourceName].dataStore;
        return Boolean(dataStore.columnIndex[targetColumnName]);
    }, columnName);
}

async function getColumnState(page: Page, chartId: string, columnName: string) {
    return await page.evaluate(({ targetChartId, targetColumnName }) => {
        const chartManager = (window as any).mdv.chartManager;
        const dataSourceName = Object.keys(chartManager.dsIndex)[0];
        const dataStore = chartManager.dsIndex[dataSourceName].dataStore;
        const tableChart = Object.values(chartManager.charts)
            .map((entry: any) => entry.chart)
            .find((chart: any) => chart.config.id === targetChartId);
        return {
            hasUnsavedChanges: chartManager.viewManager?.hasUnsavedChanges?.() ?? null,
            metadata:
                dataStore.getAllColumnsMetadata().find(
                    (column: any) => column.field === targetColumnName,
                ) ?? null,
            stillVisible: Boolean(dataStore.columnIndex[targetColumnName]),
            tableStillUsesColumn:
                Array.isArray(tableChart?.config?.param) &&
                tableChart.config.param.includes(targetColumnName),
        };
    }, { targetChartId: chartId, targetColumnName: columnName });
}

async function expectFeedbackAlert(page: Page, title: string, message: string) {
    const alert = page
        .getByRole("dialog")
        .filter({ has: page.getByText(title, { exact: true }) });
    await expect(alert).toBeVisible();
    await expect(alert.getByText(title, { exact: true })).toBeVisible();
    await expect(alert.getByText(message, { exact: true })).toBeVisible();
    return alert;
}

async function closeFeedbackAlert(page: Page) {
    const alert = page
        .getByRole("dialog")
        .filter({ has: page.getByText("Delete Column Error", { exact: true }) });
    await alert.getByRole("button", { name: "close", exact: true }).click();
    await expect(alert).not.toBeVisible();
}

async function injectSoftDeleteFailure(page: Page, columnName: string) {
    await page.evaluate((targetColumnName) => {
        const chartManager = (window as any).mdv.chartManager;
        const dataSourceName = Object.keys(chartManager.dsIndex)[0];
        const dataStore = chartManager.dsIndex[dataSourceName].dataStore;
        const originalSoftDeleteColumn = dataStore.softDeleteColumn.bind(dataStore);
        dataStore.softDeleteColumn = (field: string, dirty = false, notify = false) => {
            if (field === targetColumnName) {
                throw new Error("Injected soft delete failure");
            }
            return originalSoftDeleteColumn(field, dirty, notify);
        };
    }, columnName);
}

async function injectSaveFailure(page: Page) {
    await page.evaluate(() => {
        const viewManager = (window as any).mdv.chartManager.viewManager;
        viewManager.saveView = async () => {
            throw new Error("Injected save failure");
        };
    });
}

async function injectSavedViewAnalysisFailure(page: Page, brokenViewName: string) {
    await page.evaluate((targetViewName) => {
        const chartManager = (window as any).mdv.chartManager;
        const viewManager = chartManager.viewManager;
        const allViews = Array.isArray(viewManager.all_views) ? viewManager.all_views : [];
        if (!allViews.includes(targetViewName)) {
            viewManager.all_views = [...allViews, targetViewName];
        }
        const originalViewLoader = chartManager.viewLoader?.bind(chartManager);
        chartManager.viewLoader = async (viewName: string) => {
            if (viewName === targetViewName) {
                throw new Error("Injected view load failure");
            }
            if (!originalViewLoader) {
                return null;
            }
            return await originalViewLoader(viewName);
        };
    }, brokenViewName);
}

async function setupSoftDeleteFixture(page: Page) {
    const projectHandle = await createTemporaryProjectViaSyntheticAnndata(page, {
        synthetic: {
            profile: "minimal",
            nCells: 200,
            nGenes: 12,
            force: true,
        },
    });
    const setupInfo = await setupReactTableProject(page);
    return { projectHandle, setupInfo };
}

test.describe("Soft Delete", () => {
    test.setTimeout(180_000);

    test("blocks deleting a column when another chart in the current view still uses it", async ({ page }) => {
        const { projectHandle, setupInfo } = await setupSoftDeleteFixture(page);
        const columnName = "soft_delete_blocked";

        try {
            await addClonedEditableColumn(
                page,
                setupInfo.tableChartId,
                columnName,
                setupInfo.cloneSourceName,
            );
            await addDependentRowChart(
                page,
                setupInfo.dataSourceName,
                columnName,
                "Soft Delete Dependency",
            );

            await requestColumnDeletion(page, setupInfo.tableChartId, columnName);

            const dialog = page.getByRole("dialog", { name: `Delete Column "${columnName}"` });
            await expect(dialog).toBeVisible();
            await expect(dialog.getByText("Deletion is blocked.")).toBeVisible();
            await expect(dialog.getByText("Current View", { exact: true })).toBeVisible();
            await expect(dialog.getByText("Soft Delete Dependency")).toBeVisible();
            await expect(
                dialog.getByRole("button", { name: "Deletion Blocked", exact: true }),
            ).toBeDisabled();

            await dialog.getByRole("button", { name: "Cancel", exact: true }).click();
            await expect(dialog).not.toBeVisible();
        } finally {
            await projectHandle.cleanup();
        }
    });

    test("soft-deletes the column, saves the view, and keeps it hidden after reload", async ({ page }) => {
        const { projectHandle, setupInfo } = await setupSoftDeleteFixture(page);
        const columnName = "soft_delete_persist";

        try {
            await createSavedEditableColumn(
                page,
                setupInfo.tableChartId,
                columnName,
                setupInfo.cloneSourceName,
            );
            await requestColumnDeletion(page, setupInfo.tableChartId, columnName);

            const dialog = page.getByRole("dialog", { name: `Delete Column "${columnName}"` });
            await expect(dialog).toBeVisible();
            await expect(
                dialog.getByText(
                    "Note: This is a soft delete. The column won't be visible in MDV, but it will still appear in exported datasource files.",
                ),
            ).toBeVisible();

            await dialog.getByRole("button", { name: "Delete Column & Save", exact: true }).click();
            await expect(dialog).not.toBeVisible({ timeout: 60_000 });

            const stateAfterDelete = await getColumnState(page, setupInfo.tableChartId, columnName);

            expect(stateAfterDelete.hasUnsavedChanges).toBe(false);
            expect(stateAfterDelete.metadata).toEqual(
                expect.objectContaining({
                    field: columnName,
                    deleted: true,
                }),
            );
            expect(stateAfterDelete.stillVisible).toBe(false);
            expect(stateAfterDelete.tableStillUsesColumn).toBe(false);

            await page.reload();
            await waitForProjectReady(page);

            const stateAfterReload = await getColumnState(page, setupInfo.tableChartId, columnName);

            expect(stateAfterReload.metadata).toEqual(
                expect.objectContaining({
                    field: columnName,
                    deleted: true,
                }),
            );
            expect(stateAfterReload.stillVisible).toBe(false);
            expect(stateAfterReload.tableStillUsesColumn).toBe(false);
        } finally {
            await projectHandle.cleanup();
        }
    });

    test("blocks deleting a column when another saved view still uses it", async ({ page }) => {
        const { projectHandle, setupInfo } = await setupSoftDeleteFixture(page);
        const columnName = "soft_delete_saved_view_blocked";
        const savedViewName = `soft-delete-view-${Date.now()}`;
        const dependencyChartTitle = "Saved View Dependency";

        try {
            await createSavedEditableColumn(
                page,
                setupInfo.tableChartId,
                columnName,
                setupInfo.cloneSourceName,
            );
            await injectSavedViewImpact(
                page,
                setupInfo.dataSourceName,
                columnName,
                savedViewName,
                dependencyChartTitle,
            );

            await requestColumnDeletion(page, setupInfo.tableChartId, columnName);

            const dialog = page.getByRole("dialog", { name: `Delete Column "${columnName}"` });
            await expect(dialog).toBeVisible();
            await expect(dialog.getByText("Deletion is blocked.")).toBeVisible();
            await expect(dialog.getByText("Other Views", { exact: true })).toBeVisible();
            await expect(dialog.getByText(savedViewName, { exact: true }).first()).toBeVisible();
            await expect(dialog.getByText(dependencyChartTitle)).toBeVisible();
            await expect(
                dialog.getByRole("button", { name: "Open View", exact: true }).first(),
            ).toBeVisible();
            await expect(
                dialog.getByRole("button", { name: "Deletion Blocked", exact: true }),
            ).toBeDisabled();
            await dialog.getByRole("button", { name: "Cancel", exact: true }).click();
            await expect(dialog).not.toBeVisible();
        } finally {
            await projectHandle.cleanup();
        }
    });

    test("completely removes a newly added unsaved column instead of tombstoning it", async ({ page }) => {
        const { projectHandle, setupInfo } = await setupSoftDeleteFixture(page);
        const columnName = "soft_delete_unsaved_new_column";

        try {
            await addClonedEditableColumn(
                page,
                setupInfo.tableChartId,
                columnName,
                setupInfo.cloneSourceName,
            );
            await requestColumnDeletion(page, setupInfo.tableChartId, columnName);

            const dialog = page.getByRole("dialog", { name: `Delete Column "${columnName}"` });
            await expect(dialog).toBeVisible();
            await dialog.getByRole("button", { name: "Delete Column & Save", exact: true }).click();
            await expect(dialog).not.toBeVisible({ timeout: 60_000 });

            const stateAfterDelete = await getColumnState(page, setupInfo.tableChartId, columnName);
            expect(stateAfterDelete.hasUnsavedChanges).toBe(false);
            expect(stateAfterDelete.metadata).toBeNull();
            expect(stateAfterDelete.stillVisible).toBe(false);
            expect(stateAfterDelete.tableStillUsesColumn).toBe(false);

            await page.reload();
            await waitForProjectReady(page);

            const stateAfterReload = await getColumnState(page, setupInfo.tableChartId, columnName);
            expect(stateAfterReload.metadata).toBeNull();
            expect(stateAfterReload.stillVisible).toBe(false);
            expect(stateAfterReload.tableStillUsesColumn).toBe(false);
        } finally {
            await projectHandle.cleanup();
        }
    });

    test("shows an error and leaves the column intact when soft delete fails", async ({ page }) => {
        const { projectHandle, setupInfo } = await setupSoftDeleteFixture(page);
        const columnName = "soft_delete_delete_error";

        try {
            await createSavedEditableColumn(
                page,
                setupInfo.tableChartId,
                columnName,
                setupInfo.cloneSourceName,
            );
            await injectSoftDeleteFailure(page, columnName);
            await requestColumnDeletion(page, setupInfo.tableChartId, columnName);

            const dialog = page.getByRole("dialog", { name: `Delete Column "${columnName}"` });
            await expect(dialog).toBeVisible();
            await dialog.getByRole("button", { name: "Delete Column & Save", exact: true }).click();

            await expectFeedbackAlert(
                page,
                "Delete Column Error",
                `Column ${columnName} could not be deleted.`,
            );
            await expect(dialog).not.toBeVisible();

            const stateAfterFailure = await getColumnState(page, setupInfo.tableChartId, columnName);
            expect(stateAfterFailure.metadata).toEqual(
                expect.objectContaining({
                    field: columnName,
                }),
            );
            expect(stateAfterFailure.metadata?.deleted).not.toBe(true);
            expect(stateAfterFailure.stillVisible).toBe(true);
            expect(stateAfterFailure.tableStillUsesColumn).toBe(true);

            await closeFeedbackAlert(page);
            await page.reload();
            await waitForProjectReady(page);
        } finally {
            await projectHandle.cleanup();
        }
    });

    test("shows an error when saved-view impact analysis fails and does not start deletion", async ({ page }) => {
        const { projectHandle, setupInfo } = await setupSoftDeleteFixture(page);
        const columnName = "soft_delete_analysis_error";
        const brokenViewName = "__soft_delete_broken_view__";

        try {
            await addClonedEditableColumn(
                page,
                setupInfo.tableChartId,
                columnName,
                setupInfo.cloneSourceName,
            );
            await injectSavedViewAnalysisFailure(page, brokenViewName);
            await requestColumnDeletion(page, setupInfo.tableChartId, columnName);

            await expectFeedbackAlert(
                page,
                "Delete Column Error",
                `Failed to check column usage in saved views: ${brokenViewName}`,
            );
            await expect(
                page.getByRole("dialog", { name: `Delete Column "${columnName}"` }),
            ).not.toBeVisible();

            const stateAfterFailure = await getColumnState(page, setupInfo.tableChartId, columnName);
            expect(stateAfterFailure.metadata).toEqual(
                expect.objectContaining({
                    field: columnName,
                }),
            );
            expect(stateAfterFailure.stillVisible).toBe(true);
            expect(stateAfterFailure.tableStillUsesColumn).toBe(true);

            await closeFeedbackAlert(page);
            await page.reload();
            await waitForProjectReady(page);
        } finally {
            await projectHandle.cleanup();
        }
    });

    test("shows an error when saving fails after a local soft delete and restores on reload", async ({ page }) => {
        const { projectHandle, setupInfo } = await setupSoftDeleteFixture(page);
        const columnName = "soft_delete_save_error";

        try {
            await createSavedEditableColumn(
                page,
                setupInfo.tableChartId,
                columnName,
                setupInfo.cloneSourceName,
            );
            await injectSaveFailure(page);
            await requestColumnDeletion(page, setupInfo.tableChartId, columnName);

            const dialog = page.getByRole("dialog", { name: `Delete Column "${columnName}"` });
            await expect(dialog).toBeVisible();
            await dialog.getByRole("button", { name: "Delete Column & Save", exact: true }).click();

            await expectFeedbackAlert(
                page,
                "Delete Column Error",
                `Column ${columnName} was deleted locally, but saving the updated view failed. Reloading the project may restore this column.`,
            );
            await expect(dialog).not.toBeVisible();

            const localStateAfterFailure = await getColumnState(page, setupInfo.tableChartId, columnName);
            expect(localStateAfterFailure.metadata).toEqual(
                expect.objectContaining({
                    field: columnName,
                    deleted: true,
                }),
            );
            expect(localStateAfterFailure.stillVisible).toBe(false);
            expect(localStateAfterFailure.tableStillUsesColumn).toBe(false);

            await page.reload();
            await waitForProjectReady(page);

            const stateAfterReload = await getColumnState(page, setupInfo.tableChartId, columnName);
            expect(stateAfterReload.metadata).toEqual(
                expect.objectContaining({
                    field: columnName,
                }),
            );
            expect(stateAfterReload.metadata?.deleted).not.toBe(true);
            expect(stateAfterReload.stillVisible).toBe(true);
            expect(stateAfterReload.tableStillUsesColumn).toBe(true);
        } finally {
            await projectHandle.cleanup();
        }
    });
});
