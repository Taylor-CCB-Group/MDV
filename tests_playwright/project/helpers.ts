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
