import test, { expect, type BrowserContext, type Page } from "@playwright/test";
import { newProjectSetup, newProjectTeardown } from "../utils/testUtils";

// todo: Try to come up with a better approach to run the tests in parellel
// Running tests serially
test.describe.configure({ mode: 'serial' });

test.describe('Chart Creation', () => {
    let projectId: string;
    let context: BrowserContext;
    let page: Page;

    // Shared browser context
    test.beforeAll(async ({ browser }) => {
        context = await browser.newContext();
        page = await context.newPage();
        projectId = await newProjectSetup(page);
    });

    test.afterAll(async () => {
        await newProjectTeardown(page, projectId);
        await context.close();
    });

    // Helper function to add a chart and assert it was created
    async function addChartAndAssert(
        page: Page,
        chartName: string,
        chartType: string
    ) {
        // Add a new chart
        await page.locator('.ciview-menu-icon').first().hover();
        await expect(page.getByRole('tooltip', { name: 'Add Chart', exact: true }).first()).toBeVisible();
        await page.locator('.ciview-menu-icon').first().click();
        await expect(page.getByRole('dialog')).toBeVisible();
        await page.getByRole('combobox', { name: 'Chart Type' }).click();
        await expect.soft(page.getByText(chartName).first()).toBeVisible();
        await page.getByRole('option', { name: chartName, exact: true }).click();
        await page.getByRole('button', { name: 'Add Chart', exact: true }).click();
        
        // Wait for chartManager and verify the chart type
        await page.waitForFunction(() => !!(window as any).mdv?.chartManager);
        await page.waitForTimeout(500);
        
        // Check if the chart type exists in the config
        const hasChart = await page.evaluate((type) => {
            const chartEntries = Object.values((window as any).mdv.chartManager.charts);
            return chartEntries.some((entry: any) => entry.chart.config.type === type);
        }, chartType);
        
        expect(hasChart).toBe(true);
    }

    test('2D Scatter Plot', async () => {
        await addChartAndAssert(page, '2D Scatter Plot', 'wgl_scatter_plot');
    });

    test('3D Scatter Plot', async () => {
        await addChartAndAssert(page, '3D Scatter Plot', 'wgl_3d_scatter_plot');
    });

    test('Histogram', async () => {
        await addChartAndAssert(page, 'Histogram', 'bar_chart');
    });

    test('Row Chart', async () => {
        await addChartAndAssert(page, 'Row Chart', 'row_chart');
    });

    test('Table', async () => {
        await addChartAndAssert(page, 'Table', 'table_chart');
    });

    test('Violin Plot', async () => {
        await addChartAndAssert(page, 'Violin Plot', 'violin_plot');
    });

    test('Box Plot', async () => {
        await addChartAndAssert(page, 'Box Plot', 'box_plot');
    });

    test('Heat Map', async () => {
        await addChartAndAssert(page, 'Heat Map', 'heat_map');
    });

    test('Density Scatter Plot', async () => {
        await addChartAndAssert(page, 'Density Scatter Plot', 'density_scatter_plot');
    });

    test('Pie Chart', async () => {
        await addChartAndAssert(page, 'Pie Chart', 'ring_chart');
    });

    test('Sankey Diagram', async () => {
        await addChartAndAssert(page, 'Sankey Diagram', 'sankey_chart');
    });

    test('Stacked Row Chart', async () => {
        await addChartAndAssert(page, 'Stacked Row Chart', 'stacked_row_chart');
    });

    test('Multi Line Chart', async () => {
        await addChartAndAssert(page, 'Multi Line Chart', 'multi_line_chart');
    });

    test('Word Cloud', async () => {
        await addChartAndAssert(page, 'Word Cloud', 'row_chart');
    });

    test('Dot Plot', async () => {
        await addChartAndAssert(page, 'Dot Plot', 'dot_plot');
    });

    test('Selection Dialog', async () => {
        await addChartAndAssert(page, 'Selection Dialog', 'selection_dialog');
    });

    test('Text Box', async () => {
        await addChartAndAssert(page, 'Text Box', 'text_box_chart');
    });

    test('Row Summary Box', async () => {
        await addChartAndAssert(page, 'Row Summary Box', 'row_summary_box');
    });

    test('Abundance Box Plot', async () => {
        await addChartAndAssert(page, 'Abundance Box Plot', 'custom_box_plot');
    });

});