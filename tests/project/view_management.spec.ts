import test, { expect, type BrowserContext, type Page } from "@playwright/test";
import { newProjectSetup, newProjectTeardown } from "../utils/testUtils";

test.describe.configure({ mode: 'serial' });

test.describe('View Management', () => {
    let projectId: string;
    let context: BrowserContext;
    let page: Page;

    test.beforeAll(async ({ browser }) => {
        context = await browser.newContext();
        page = await context.newPage();
        projectId = await newProjectSetup(page);
    });

    test.afterAll(async () => {
        await newProjectTeardown(page, projectId);
        await context.close();
    });

    // Helper function to get current view
    async function getCurrentView(page: Page): Promise<string> {
        await page.waitForFunction(() => !!(window as any).mdv?.chartManager?.viewManager);
        return await page.evaluate(() => (window as any).mdv.chartManager.viewManager.current_view);
    }

    // Helper function to get all views
    async function getAllViews(page: Page): Promise<string[]> {
        await page.waitForFunction(() => !!(window as any).mdv?.chartManager?.viewManager);
        return await page.evaluate(() => (window as any).mdv.chartManager.viewManager.all_views);
    }

    // Helper function to create a new view
    async function createView(page: Page, viewName: string) {
        await page.getByRole('button', { name: /create new view/i }).click();
        await expect(page.getByRole('dialog')).toBeVisible();
        await page.getByRole('textbox', { name: 'Name', exact: true }).fill(viewName);
        await page.getByRole('button', { name: /create view/i, exact: true }).click();
        await expect(page.getByRole('dialog')).not.toBeVisible();
        await page.waitForTimeout(1000);
    }

    test('verify default view is loaded', async () => {
        const currentView = await getCurrentView(page);
        expect(currentView).toBeDefined();
        expect(currentView).toBe('default');
    });

    test('create new view', async () => {
        const initialViews = await getAllViews(page);
        const initialCount = initialViews.length;

        // Create new view
        const newViewName = `test_view_${Date.now()}`;
        await createView(page, newViewName);

        // Verify view was created
        const currentViews = await getAllViews(page);
        expect(currentViews.length).toBe(initialCount + 1);
        expect(currentViews).toContain(newViewName);

        // Verify current view changed to new view
        const currentView = await getCurrentView(page);
        expect(currentView).toBe(newViewName);
    });

    test('switch view using view selector dropdown', async () => {
        const currentView = await getCurrentView(page);

        const newViewName = `switch_test_view_${Date.now()}`;
        await createView(page, newViewName);

        // Switch back to original view
        await page.getByLabel(/select view/i).click();
        await page.getByRole('option', { name: currentView }).click();
        await page.waitForTimeout(1000);

        // Switch to target view
        await page.getByLabel(/select view/i).click();
        await page.getByRole('option', { name: newViewName }).click();
        await page.waitForTimeout(1000);

        // Verify view changed
        const newCurrentView = await getCurrentView(page);
        expect(newCurrentView).toBe(newViewName);
    });

    test('access view from view gallery', async () => {
        const newViewName = `gallery_test_view_${Date.now()}`;
        await createView(page, newViewName);
        await page.getByRole('button', { name: 'Browse View Gallery' }).click();
        await expect(page.getByRole('dialog')).toBeVisible();

        await page.getByText(newViewName).click();

        await page.waitForTimeout(1000);

        // Verify view changed
        const newCurrentView = await getCurrentView(page);
        expect(newCurrentView).toBe(newViewName);
    });

    test('save current view', async () => {
        const initialCharts = (await page.locator('.ciview-chart-panel').all()).length;
        // Add a chart to make a change
        await page.locator('.ciview-menu-icon').first().hover();
        await expect(page.getByRole('tooltip', { name: 'Add Chart', exact: true }).first()).toBeVisible();
        await page.locator('.ciview-menu-icon').first().click();
        await expect(page.getByRole('dialog')).toBeVisible();
        await page.getByRole('combobox', { name: 'Chart Type' }).click();
        await page.getByRole('option', { name: '2D Scatter Plot', exact: true }).click();
        await page.getByRole('button', { name: 'Add Chart', exact: true }).click();

        await page.waitForTimeout(1000);

        // Click save button
        await page.getByRole('button', { name: 'Save View', exact: true }).click();

        // Wait for save to complete
        await page.waitForTimeout(1000);

        const currentCharts = (await page.locator('.ciview-chart-panel').all()).length;

        // Check if the new chart exists
        expect(currentCharts - initialCharts).toBe(1);
    });

    test('save view as new name', async () => {
        const beforeViews = await getAllViews(page);

        // Click save as button
        await page.getByRole('button', { name: 'Save View As' }).click();
        await expect(page.getByRole('dialog')).toBeVisible();

        // Enter view name
        const newViewName = `saved_as_view_${Date.now()}`;
        await page.getByLabel(/enter new view name/i).fill(newViewName);
        await page.getByRole('button', { name: 'Save As', exact: true }).click();

        await expect(page.getByRole('dialog')).not.toBeVisible();
        await page.waitForTimeout(1000);

        // Verify new view was created
        const afterViews = await getAllViews(page);
        expect(afterViews.length).toBe(beforeViews.length + 1);
        expect(afterViews).toContain(newViewName);

        // Verify current view changed to new view
        const newCurrentView = await getCurrentView(page);
        expect(newCurrentView).toBe(newViewName);
    });

    test('delete a view', async () => {
        // Create a view to delete
        const newViewName = `delete_test_view_${Date.now()}`;
        await createView(page, newViewName);

        const currentViews = await getAllViews(page);

        const currentView = await getCurrentView(page);

        // Change to new view if not already the current view
        if (currentView !== newViewName) {
            await page.getByLabel(/select view/i).click();
            await page.getByRole('option', { name: newViewName }).click();
            await page.waitForTimeout(1000);
        }

        await page.getByRole('button', { name: /delete current view/i }).click();

        await expect(page.getByRole('dialog')).toBeVisible();
        await page.getByRole('button', { name: 'Yes', exact: true }).click();

        await expect(page.getByRole('dialog')).not.toBeVisible();
        await page.waitForTimeout(1000);

        // Verify view was deleted
        const afterViews = await getAllViews(page);
        expect(afterViews.length).toBe(currentViews.length - 1);
        expect(afterViews).not.toContain(newViewName);

        // Verify current view changed to a different view
        const newCurrentView = await getCurrentView(page);
        expect(newCurrentView).not.toBe(newViewName);
        expect(afterViews).toContain(newCurrentView);
    });

    test('view persists after page reload', async () => {
        const currentView = await getCurrentView(page);

        // Reload the page
        await page.reload();
        await page.waitForLoadState('networkidle');
        await page.waitForFunction(() => !!(window as any).mdv?.chartManager);

        // Verify view is still the same
        const viewAfterReload = await getCurrentView(page);
        expect(viewAfterReload).toBe(currentView);
    });
});

