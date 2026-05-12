import test, { expect, type Page } from "@playwright/test";
import {
    createTemporaryProjectViaSyntheticAnndata,
    type SyntheticAnndataTemporaryProjectHandle,
    waitForProjectReady,
} from "../../utils/projectFixtures";
import { saveCurrentView } from "../../utils/helpers";

const TABLE_CHART_TITLE = "Rename Column Table";

type ProjectSetupInfo = {
    cloneSourceName: string;
    dataSourceName: string;
    tableChartId: string;
};

/** Same `GET …/datasources.json` URL the debug toolbar uses (see `DebugButton.tsx`). */
function datasourcesJsonUrlFromProjectPage(pageUrl: string): string {
    const u = new URL(pageUrl);
    const root = `${u.origin}${u.pathname}`.replace(/\/$/, "");
    return `${root}/datasources.json`;
}

type DatasourceJsonEntry = {
    name?: string;
    columns?: Array<{ field?: string; name?: string }>;
};

async function getColumnNameFromDatasourcesJson(
    page: Page,
    dataSourceName: string,
    columnField: string,
): Promise<string | undefined> {
    const url = datasourcesJsonUrlFromProjectPage(page.url());
    const response = await page.request.get(url);
    expect(response.ok(), `GET ${url} → ${response.status()}`).toBeTruthy();
    const datasources = (await response.json()) as DatasourceJsonEntry[];
    const ds =
        datasources.find((entry) => entry.name === dataSourceName) ??
        (datasources.length > 0 ? datasources[0] : undefined);
    const col = ds?.columns?.find((c) => c.field === columnField);
    return col?.name !== undefined && col?.name !== null ? String(col.name) : undefined;
}

/** Runtime display name (not written to `datasources.json` until save). */
async function getColumnDisplayName(page: Page, columnField: string): Promise<string> {
    return await page.evaluate((field) => {
        const chartManager = (window as any).mdv.chartManager;
        const dataSourceName = Object.keys(chartManager.dsIndex)[0];
        const dataStore = chartManager.dsIndex[dataSourceName].dataStore;
        return String(dataStore.columnIndex[field]?.name ?? "");
    }, columnField);
}

async function setupReactTableProject(page: Page): Promise<ProjectSetupInfo> {
    return await page.evaluate(async (tableTitle) => {
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
        const tableChartId = `rename-col-table-${Date.now()}`;

        await chartManager.addChart(
            dataSourceName,
            {
                id: tableChartId,
                title: tableTitle,
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
    }, TABLE_CHART_TITLE);
}

async function addClonedEditableColumn(
    page: Page,
    chartId: string,
    columnName: string,
    cloneSourceName: string,
    tableTitle: string,
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

    await page.waitForFunction(
        ({ targetColumnName, chartTitle }) => {
            const chartManager = (window as any).mdv.chartManager;
            const dataSourceName = Object.keys(chartManager.dsIndex)[0];
            const dataStore = chartManager.dsIndex[dataSourceName].dataStore;
            const tableChart = Object.values(chartManager.charts)
                .map((entry: any) => entry.chart)
                .find((chart: any) => chart.config.title === chartTitle);
            return Boolean(
                dataStore.columnIndex[targetColumnName] &&
                    Array.isArray(tableChart?.config?.param) &&
                    tableChart.config.param.includes(targetColumnName),
            );
        },
        { targetColumnName: columnName, chartTitle: tableTitle },
    );
}

async function requestRenameColumnDialog(page: Page, chartId: string, columnField: string) {
    await page.evaluate(({ targetChartId, targetColumnField }) => {
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
            slickGrid.getColumns().find((entry: any) => entry.field === targetColumnField) ?? {
                field: targetColumnField,
                id: targetColumnField,
            };
        pubSub.publish("onHeaderMenuCommand", {
            column,
            command: "rename-column",
        });
    }, { targetChartId: chartId, targetColumnField: columnField });
}

async function submitRenameColumnDialog(page: Page, initialColumnName: string, newDisplayName: string) {
    const renameDialog = page.getByRole("dialog", { name: `Rename "${initialColumnName}"` });
    await expect(renameDialog).toBeVisible();
    await renameDialog.getByLabel("Column Name").fill(newDisplayName);
    await renameDialog.getByRole("button", { name: "Rename", exact: true }).click();
    await expect(renameDialog).not.toBeVisible();
}

async function setupRenameFixture(page: Page): Promise<{
    projectHandle: SyntheticAnndataTemporaryProjectHandle;
    setupInfo: ProjectSetupInfo;
}> {
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

test.describe("Rename table column", () => {
    test.setTimeout(180_000);

    test("persists a renamed newly added column after save and reload", async ({ page }) => {
        const { projectHandle, setupInfo } = await setupRenameFixture(page);
        const columnField = `rename_persist_${Date.now()}`;
        const newDisplayName = `Renamed ${Date.now()}`;

        try {
            await addClonedEditableColumn(
                page,
                setupInfo.tableChartId,
                columnField,
                setupInfo.cloneSourceName,
                TABLE_CHART_TITLE,
            );
            expect(await getColumnDisplayName(page, columnField)).toBe(columnField);

            await requestRenameColumnDialog(page, setupInfo.tableChartId, columnField);
            await submitRenameColumnDialog(page, columnField, newDisplayName);
            expect(await getColumnDisplayName(page, columnField)).toBe(newDisplayName);

            await saveCurrentView(page);
            await expect
                .poll(async () =>
                    getColumnNameFromDatasourcesJson(page, setupInfo.dataSourceName, columnField),
                )
                .toBe(newDisplayName);

            await page.reload();
            await waitForProjectReady(page);

            expect(
                await getColumnNameFromDatasourcesJson(page, setupInfo.dataSourceName, columnField),
            ).toBe(newDisplayName);
        } finally {
            await projectHandle.cleanup();
        }
    });

    test("drops an unsaved rename after reload when the column was saved earlier", async ({ page }) => {
        const { projectHandle, setupInfo } = await setupRenameFixture(page);
        const columnField = `rename_revert_${Date.now()}`;
        const unsavedRename = `Should Not Persist ${Date.now()}`;

        try {
            await addClonedEditableColumn(
                page,
                setupInfo.tableChartId,
                columnField,
                setupInfo.cloneSourceName,
                TABLE_CHART_TITLE,
            );
            await saveCurrentView(page);
            await expect
                .poll(async () =>
                    getColumnNameFromDatasourcesJson(page, setupInfo.dataSourceName, columnField),
                )
                .toBe(columnField);

            await page.reload();
            await waitForProjectReady(page);

            expect(
                await getColumnNameFromDatasourcesJson(page, setupInfo.dataSourceName, columnField),
            ).toBe(columnField);

            await requestRenameColumnDialog(page, setupInfo.tableChartId, columnField);
            await submitRenameColumnDialog(page, columnField, unsavedRename);
            expect(await getColumnDisplayName(page, columnField)).toBe(unsavedRename);
            expect(
                await getColumnNameFromDatasourcesJson(page, setupInfo.dataSourceName, columnField),
            ).toBe(columnField);

            await page.reload();
            await waitForProjectReady(page);

            expect(
                await getColumnNameFromDatasourcesJson(page, setupInfo.dataSourceName, columnField),
            ).toBe(columnField);
        } finally {
            await projectHandle.cleanup();
        }
    });
});
