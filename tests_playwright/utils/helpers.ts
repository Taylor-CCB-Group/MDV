import { expect, type Page } from "@playwright/test";

export const GENERATED_MOCK_PROJECT = {
    nCells: 8,
    nGenes: 12,
} as const;

export type ChartSummary = {
    id: string;
    title: string;
    type: string;
};

export type AvailableChartType = {
    chartName: string;
    expectedType: string;
};

type ChartEntry = {
    chart?: {
        config?: {
            id?: unknown;
            title?: unknown;
            type?: unknown;
        };
    };
};

type RawChartConfig = {
    id?: unknown;
    title?: unknown;
    type?: unknown;
};

type ChartConfig = {
    id: unknown;
    title?: unknown;
    type?: unknown;
};

export async function getChartSummaries(page: Page): Promise<ChartSummary[]> {
    return await page.evaluate(() => {
        const chartEntries = Object.values((window as any).mdv.chartManager.charts);
        const isChartEntry = (value: unknown): value is ChartEntry => {
            if (typeof value !== "object" || value === null) {
                return false;
            }
            return "chart" in value;
        };
        const hasId = (config: RawChartConfig | undefined): config is ChartConfig =>
            config !== undefined && config.id !== undefined && config.id !== null;
        return chartEntries
            .filter(isChartEntry)
            .map((entry) => entry.chart?.config)
            .filter(hasId)
            .map((config) => ({
                id: String(config.id),
                title: String(config.title ?? ""),
                type: String(config.type ?? ""),
            }));
    });
}

export async function getAvailableAddableChartTypes(page: Page): Promise<AvailableChartType[]> {
    return await page.evaluate(() => {
        const chartTypes = (window as any).mdv?.chartTypes ?? {};
        const chartManager = (window as any).mdv?.chartManager;
        const dataSourceName = Object.keys(chartManager?.dsIndex ?? {})[0];
        const dataStore = dataSourceName ? chartManager.dsIndex[dataSourceName]?.dataStore : undefined;
        if (!dataStore) {
            return [];
        }

        const results: Array<{ chartName: string; expectedType: string }> = [];
        const seenNames = new Set<string>();

        for (const [typeKey, chartType] of Object.entries(chartTypes)) {
            if (!chartType || typeof chartType !== "object") {
                continue;
            }
            const candidate = chartType as {
                allow_user_add?: boolean;
                name?: unknown;
                required?: string[] | ((ds: unknown) => boolean);
            };
            if (candidate.allow_user_add === false) {
                continue;
            }

            const required = candidate.required;
            if (required) {
                if (typeof required === "function") {
                    if (!required(dataStore)) {
                        continue;
                    }
                } else if (!required.every((key) => Boolean((dataStore as Record<string, unknown>)[key]))) {
                    continue;
                }
            }

            const chartName = String(candidate.name ?? "");
            if (!chartName || seenNames.has(chartName)) {
                continue;
            }
            seenNames.add(chartName);
            results.push({
                chartName,
                expectedType: String(typeKey),
            });
        }

        results.sort((a, b) => a.chartName.localeCompare(b.chartName));
        return results;
    });
}

export async function addChartViaUi(page: Page, chartName: string, title: string) {
    await page.locator(".ciview-menu-icon").first().hover();
    await expect(page.getByRole("tooltip", { name: "Add Chart", exact: true }).first()).toBeVisible();
    await page.locator(".ciview-menu-icon").first().click();

    const dialog = page.getByRole("dialog");
    await expect(dialog).toBeVisible();

    await dialog.getByRole("combobox", { name: "Chart Type" }).click();
    await page.getByRole("option", { name: chartName, exact: true }).click();
    await dialog.getByRole("textbox", { name: "Title", exact: true }).fill(title);
    await dialog.getByRole("button", { name: "Add Chart", exact: true }).click();

    await expect(dialog).not.toBeVisible();
}

export function getChartPanelByTitle(page: Page, title: string) {
    return page
        .locator(".ciview-chart-panel")
        .filter({ has: page.locator(".ciview-chart-title", { hasText: title }) })
        .first();
}

export async function expectChartPanelHasNoError(page: Page, title: string) {
    const panel = getChartPanelByTitle(page, title);
    await expect(panel).toBeVisible();
    await expect(panel.getByText(/Error Occurred|An error occurred while creating the chart/i)).toHaveCount(
        0,
    );
}

export async function expectChartCreationErrorVisible(page: Page, message?: RegExp | string) {
    await expect(page.getByRole("button", { name: "ERROR: Click to view details", exact: true })).toBeVisible();
    await page.getByRole("button", { name: "ERROR: Click to view details", exact: true }).click();
    if (message) {
        await expect(page.getByText(message).first()).toBeVisible();
    }
}

export async function openChartSettingsDialog(page: Page, title: string) {
    const panel = getChartPanelByTitle(page, title);
    await expect(panel).toBeVisible();
    await panel.locator(".ciview-chart-menuspace > span").nth(2).click();
    await expect(page.locator(".ciview-dlg-header-text", { hasText: `Settings: ${title}` })).toBeVisible();
    await expect(page.getByLabel("Search Settings by Folder or Name")).toBeVisible();
}

export async function openChartDebugDialog(page: Page, title: string) {
    const panel = getChartPanelByTitle(page, title);
    await expect(panel).toBeVisible();
    await panel.locator(".ciview-chart-menuspace > span").nth(1).click();
    await page.getByText("debug / report chart", { exact: true }).click();
    await expect(page.locator(".ciview-dlg-header-text", { hasText: `Debug / Report ${title}` })).toBeVisible();
}

export async function closeTopDialog(page: Page) {
    await page.locator(".ciview-dlg-close-icon").last().click();
}

export async function waitForChartByTitle(page: Page, title: string) {
    await expect
        .poll(async () => {
            const charts = await getChartSummaries(page);
            return charts.find((chart) => chart.title === title) ?? null;
        })
        .not.toBeNull();
}

export async function waitForViewUnsavedState(page: Page, expected: boolean) {
    await page.waitForFunction((isUnsaved) => {
        const viewManager = (window as any).mdv?.chartManager?.viewManager;
        return viewManager?.hasUnsavedChanges?.() === isUnsaved;
    }, expected);
}

export async function saveCurrentView(page: Page) {
    await page.getByRole("button", { name: "Save View", exact: true }).click();
    await waitForViewUnsavedState(page, false);
}

export async function getCurrentView(page: Page): Promise<string> {
    await page.waitForFunction(() => Boolean((window as any).mdv?.chartManager?.viewManager));
    return await page.evaluate(() => String((window as any).mdv.chartManager.viewManager.current_view));
}

export async function getAllViews(page: Page): Promise<string[]> {
    await page.waitForFunction(() => Boolean((window as any).mdv?.chartManager?.viewManager));
    return await page.evaluate(() => {
        const allViews = (window as any).mdv.chartManager.viewManager.all_views;
        return Array.isArray(allViews) ? allViews.map((view) => String(view)) : [];
    });
}

export async function createViewViaUi(page: Page, viewName: string) {
    await page.getByRole("button", { name: /create new view/i }).click();
    const dialog = page.getByRole("dialog");
    await expect(dialog).toBeVisible();
    await dialog.getByRole("textbox", { name: "Name", exact: true }).fill(viewName);
    await dialog.getByRole("button", { name: /create view/i, exact: true }).click();
    await expect(dialog).not.toBeVisible();
    await expect.poll(async () => await getCurrentView(page)).toBe(viewName);
}

export async function selectViewViaUi(page: Page, viewName: string) {
    await page.getByRole("combobox", { name: "Select View" }).click();
    await page.getByRole("option", { name: viewName, exact: true }).click();
    await expect.poll(async () => await getCurrentView(page)).toBe(viewName);
}

export type ScatterDensityMode = "grid" | "overlay";

export async function patchScatterChartConfig(
    page: Page,
    title: string,
    patch: { density_mode?: ScatterDensityMode; densityFields?: string[]; type?: string },
) {
    await page.evaluate(
        ({ chartTitle, configPatch }) => {
            const chartEntries = Object.values((window as any).mdv.chartManager.charts) as Array<{
                chart?: { config?: Record<string, unknown> };
            }>;
            const entry = chartEntries.find((item) => item.chart?.config?.title === chartTitle);
            const config = entry?.chart?.config;
            if (!config) {
                throw new Error(`Chart not found: ${chartTitle}`);
            }
            if (configPatch.type) {
                config.type = configPatch.type;
            }
            if (configPatch.density_mode) {
                config.density_mode = configPatch.density_mode;
                if (configPatch.density_mode === "grid") {
                    config.contour_fill = true;
                }
            }
            if (configPatch.densityFields) {
                config.densityFields = configPatch.densityFields;
            }
        },
        { chartTitle: title, configPatch: patch },
    );
}

export async function getNumericColumnNamesForChart(page: Page, title: string, limit = 2): Promise<string[]> {
    return await page.evaluate(
        ({ chartTitle, maxCount }) => {
            const chartEntries = Object.values((window as any).mdv.chartManager.charts) as Array<{
                chart?: {
                    config?: { title?: string };
                    dataStore?: { columnIndex?: Record<string, { datatype?: string }> };
                };
            }>;
            const entry = chartEntries.find((item) => item.chart?.config?.title === chartTitle);
            const columnIndex = entry?.chart?.dataStore?.columnIndex;
            if (!columnIndex) {
                return [];
            }
            return Object.entries(columnIndex)
                .filter(([, column]) => column.datatype === "double" || column.datatype === "integer")
                .map(([name]) => name)
                .slice(0, maxCount);
        },
        { chartTitle: title, maxCount: limit },
    );
}

export function getSelectionToolbar(page: Page, chartTitle: string) {
    const panel = getChartPanelByTitle(page, chartTitle);
    return panel.getByRole("group", { name: "choose tool for manipulating view or selection" });
}

export type RowsAsColumnsSubgroupInfo = {
    key: string;
    label: string;
};

type PlaywrightDsIndexEntry = {
    dataStore?: {
        links?: Record<
            string,
            {
                rows_as_columns?: { subgroups?: Record<string, { label?: string }> };
            }
        >;
    };
};

/** Subgroup keys on the cells→genes rows_as_columns link (requires extra-expression-layers project). */
export async function getRowsAsColumnsSubgroupKeys(page: Page): Promise<RowsAsColumnsSubgroupInfo[]> {
    return await page.evaluate(() => {
        const chartManager = (window as any).mdv?.chartManager;
        const dsIndex = (chartManager?.dsIndex ?? {}) as Record<string, PlaywrightDsIndexEntry>;
        const cellsStore =
            chartManager?.dsIndex?.cells?.dataStore ??
            Object.values(dsIndex).find((entry) => entry?.dataStore?.links)?.dataStore;
        if (!cellsStore?.links) {
            return [];
        }
        for (const linkedDsName of Object.keys(cellsStore.links)) {
            const subgroups = cellsStore.links[linkedDsName]?.rows_as_columns?.subgroups;
            if (!subgroups) {
                continue;
            }
            return Object.entries(subgroups).map(([key, info]) => ({
                key,
                label: String((info as { label?: string })?.label ?? key),
            }));
        }
        return [];
    });
}

export async function getFirstLinkedGeneLabel(page: Page): Promise<string> {
    return await page.evaluate(async () => {
        const chartManager = (window as any).mdv?.chartManager;
        const cellsStore =
            chartManager?.dsIndex?.cells?.dataStore ??
            chartManager.dsIndex[Object.keys(chartManager.dsIndex)[0]]?.dataStore;
        if (!cellsStore?.links) {
            throw new Error("No rows_as_columns links on cells datasource");
        }
        for (const linkedDsName of Object.keys(cellsStore.links)) {
            const rac = cellsStore.links[linkedDsName]?.rows_as_columns as { name_column?: string } | undefined;
            if (!rac?.name_column) {
                continue;
            }
            const linkedDs = chartManager.dataSources.find((ds: { name: string }) => ds.name === linkedDsName);
            const nameColumn = linkedDs?.dataStore?.columnIndex?.[rac.name_column];
            if (!nameColumn) {
                continue;
            }
            if (!nameColumn.values && linkedDs?.dataStore?.loadColumn) {
                await linkedDs.dataStore.loadColumn(nameColumn.field);
            }
            const values = nameColumn.values ?? linkedDs?.dataStore?.getColumnValues?.(rac.name_column);
            if (Array.isArray(values) && values.length > 0) {
                return String(values[0]);
            }
        }
        throw new Error("No linked gene labels available");
    });
}

export async function addDotPlotChartViaApi(page: Page, title: string): Promise<string> {
    return await page.evaluate(async (chartTitle) => {
        const chartManager = (window as any).mdv.chartManager;
        const dataSourceName = chartManager.dsIndex.cells ? "cells" : Object.keys(chartManager.dsIndex)[0];
        const dataStore = chartManager.dsIndex[dataSourceName].dataStore;
        const categoryColumn =
            dataStore.columns.find(
                (column: { datatype?: string; subgroup?: unknown; field?: string }) =>
                    column.datatype === "text" && !column.subgroup && column.field !== "__index__",
            )?.field ?? "cell_type";
        const chartId = `pw-dot-plot-${Date.now()}`;
        await chartManager.addChart(
            dataSourceName,
            {
                id: chartId,
                title: chartTitle,
                legend: "",
                type: "dot_plot",
                param: [categoryColumn],
                size: [700, 420],
                position: [20, 20],
            },
            true,
        );
        return chartId;
    }, title);
}

export async function searchChartSettings(page: Page, query: string) {
    const search = page.getByLabel("Search Settings by Folder or Name");
    await expect(search).toBeVisible();
    await search.fill(query);
    await expect(page.getByText(query, { exact: true }).first()).toBeVisible();
}

function getSettingRow(page: Page, paramLabel: string) {
    return page
        .locator(".grid")
        .filter({ has: page.getByText(paramLabel, { exact: true }) })
        .first();
}

export function getLinkColumnSettingRow(page: Page, paramLabel: string) {
    const row = getSettingRow(page, paramLabel);
    return { row, subgroupTablist: row.getByRole("tablist").nth(1) };
}

/** Open the column-selection "link" tab and pick one gene from a rows_as_columns subgroup. */
export async function selectLinkedGeneInChartSettings(
    page: Page,
    options: {
        paramLabel: string;
        subgroupTabLabel: string;
        geneLabel: string;
    },
) {
    await searchChartSettings(page, options.paramLabel);
    const { row, subgroupTablist } = getLinkColumnSettingRow(page, options.paramLabel);
    await expect(row).toBeVisible();
    await row.getByRole("tab", { name: "link", exact: true }).click();
    await subgroupTablist.getByRole("tab", { name: options.subgroupTabLabel, exact: true }).click();

    const combobox = row.getByRole("combobox").first();
    await combobox.click();
    await combobox.fill(options.geneLabel);
    const listbox = page.getByRole("listbox");
    await expect(listbox).toBeVisible();
    const option = listbox
        .getByRole("option")
        .filter({ hasText: options.geneLabel })
        .first();
    await expect(option).toBeVisible();
    await option.click();
}

/** True when Playwright targets a worktree Vite dev server (subgroup tab UI is not on the 5055 image build). */
export function isWorktreeFrontendForPlaywright(): boolean {
    const baseUrl = process.env.TEST_BASE_URL ?? "http://localhost:5055/";
    return /5170|5173/.test(baseUrl);
}

export async function setDotPlotXAxisLinkField(
    page: Page,
    chartTitle: string,
    subgroupKey: string,
    geneLabel: string,
    geneIndex: number,
) {
    await page.evaluate(
        ({ chartTitle: title, subgroupKey: sg, geneLabel: label, geneIndex: index }) => {
            const fieldName = `${sg}|${label} (${sg})|${index}`;
            const chartEntries = Object.values((window as any).mdv.chartManager.charts) as Array<{
                chart?: {
                    config?: { title?: string; param?: unknown[] };
                    _setFields?: (fields: string[]) => void;
                };
            }>;
            const entry = chartEntries.find((item) => item.chart?.config?.title === title);
            const chart = entry?.chart;
            const config = chart?.config;
            const category = config?.param?.[0];
            if (!chart || !config || typeof category !== "string") {
                throw new Error(`Dot plot chart not found or missing category: ${title}`);
            }
            config.param = [category, fieldName];
        },
        { chartTitle, subgroupKey, geneLabel, geneIndex },
    );
    await expect
        .poll(async () => {
            const fields = await getDotPlotXAxisFields(page, chartTitle);
            return (fields[0] ?? "").startsWith(`{subgroupKey}|`);
        })
        .toBe(true);
    await waitForViewUnsavedState(page, true);
}

export async function getDotPlotXAxisFields(page: Page, chartTitle: string): Promise<string[]> {
    return await page.evaluate((title) => {
        const chartEntries = Object.values((window as any).mdv.chartManager.charts) as Array<{
            chart?: { config?: { title?: string; param?: unknown[] } };
        }>;
        const entry = chartEntries.find((item) => item.chart?.config?.title === title);
        const param = entry?.chart?.config?.param;
        if (!Array.isArray(param)) {
            return [];
        }
        return param.slice(1).filter((value): value is string => typeof value === "string");
    }, chartTitle);
}
