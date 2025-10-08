import { expect, type Page } from '@playwright/test';
import { gotoPath, mockApiRoot } from './routes';
import path from 'node:path';

export const newProjectSetup = async (page: Page) => {
    await mockApiRoot(page);
    await gotoPath(page);
    await expect(page.getByTestId(/project_card/i).first()).toBeVisible();

    // Create new project
    await page.getByText('Create new project').click();
    await page.waitForURL('**/project/**');

    // Extract project id
    const url = new URL(page.url());
    const dir = url.searchParams.get('dir'); // Gets '/project/375' for ?dir=/project/375
    const newProjectId = dir
        ? dir.split('/').pop() // Extract from query param: ?dir=/project/375
        : url.pathname.split('/').pop(); // Extract from pathname: /project/375
    // Upload dataset
    await page.waitForTimeout(1000);
    await page.locator('input[type="file"]').setInputFiles(path.join(__dirname, '..', 'test-data', 'scanpy-pbmc3k.h5ad'));
    await page.getByRole('button', { name: /upload/i }).click();
    await page.getByRole('button', { name: /refresh page/i }).click();
    await expect(page.locator('.ciview-contentDiv').first()).toBeVisible();

    // Wait for charts to load
    await page.waitForTimeout(3000);
    return newProjectId as string;
};

export const newProjectTeardown = async (page: Page, projectId: string) => {
    await mockApiRoot(page);
    await gotoPath(page);
    await expect(page.getByTestId(/project_card/i).first()).toBeVisible();

    // Delete the new project
    const beforeCount = (await page.getByTestId(/project_card/i).all()).length;
    const deleteResponse = await page.request.delete(`delete_project/${projectId}`);
    expect(deleteResponse.ok()).toBeTruthy();

    // Refresh the page so the grid reflects the deletion
    await page.reload();
    await page.waitForLoadState('networkidle');
    const afterCount = (await page.getByTestId(/project_card/i).all()).length;
    expect(beforeCount - afterCount).toBe(1);
};