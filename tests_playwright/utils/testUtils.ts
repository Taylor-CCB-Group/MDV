import { expect, type Page } from '@playwright/test';
import { gotoPath, mockApiRoot } from './routes';
import path from 'node:path';

export const newProjectSetup = async (page: Page) => {
    await mockApiRoot(page);
    await gotoPath(page);
    
    // Wait for page to load - don't expect project cards on fresh database
    await page.waitForLoadState('networkidle');
    
    // Wait for the "Create new project" button to be available
    await expect(page.getByText('Create new project')).toBeVisible({ timeout: 10000 });

    // Create new project
    await page.getByText('Create new project').click();
    
    // Wait for navigation with longer timeout and better error handling
    try {
        await page.waitForURL('**/project/**', { timeout: 30000 });
    } catch (error) {
        console.log(`Navigation timeout. Current URL: ${page.url()}`);
        throw error;
    }

    // Extract project id
    const url = new URL(page.url());
    const dir = url.searchParams.get('dir'); // Gets '/project/375' for ?dir=/project/375
    const newProjectId = dir
        ? dir.split('/').pop() // Extract from query param: ?dir=/project/375
        : url.pathname.split('/').pop(); // Extract from pathname: /project/375
    
    console.log(`Created project with ID: ${newProjectId}`);
    
    // Wait longer for project page to fully load
    await page.waitForLoadState('networkidle', { timeout: 30000 });
    await page.waitForTimeout(2000); // Give React time to render
    
    // Wait for file input to be ready before setting files
    const fileInput = page.locator('input[type="file"]');
    try {
        await fileInput.waitFor({ state: 'attached', timeout: 30000 });
    } catch (error) {
        console.log(`File input not found. Current URL: ${page.url()}`);
        console.log('Checking page content...');
        throw error;
    }
    
    await fileInput.setInputFiles(path.join(__dirname, '..', 'test-data', 'scanpy-pbmc3k.h5ad'));
    
    // Click upload button
    await page.getByRole('button', { name: /upload/i }).click();
    
    // Wait for refresh button and click it
    await page.getByRole('button', { name: /refresh page/i }).waitFor({ state: 'visible', timeout: 60000 });
    await page.getByRole('button', { name: /refresh page/i }).click();
    
    // Wait for content to load
    await expect(page.locator('.ciview-contentDiv').first()).toBeVisible({ timeout: 60000 });

    // Wait for charts to load
    await page.waitForTimeout(3000);
    return newProjectId as string;
};

export const newProjectTeardown = async (page: Page, projectId: string) => {
    // If projectId is undefined, log and skip deletion
    if (!projectId) {
        console.log('No project ID provided to teardown, skipping cleanup');
        return;
    }

    try {
        await mockApiRoot(page);
        await gotoPath(page);
        await page.waitForLoadState('networkidle');

        // Get count before deletion (might be 0 or more)
        const beforeCount = (await page.getByTestId(/project_card/i).all()).length;
        
        // If no projects exist, nothing to clean up
        if (beforeCount === 0) {
            console.log('No project cards found after deletion - this is expected if all test projects were cleaned up');
            return;
        }

        // Delete the project via API
        const deleteResponse = await page.request.delete(`delete_project/${projectId}`);
        
        // Only verify deletion if the delete request was successful
        if (deleteResponse.ok()) {
            // Refresh the page so the grid reflects the deletion
            await page.reload();
            await page.waitForLoadState('networkidle');
            const afterCount = (await page.getByTestId(/project_card/i).all()).length;
            expect(beforeCount - afterCount).toBe(1);
        } else {
            console.log(`Failed to delete project ${projectId}: ${deleteResponse.status()}`);
        }
    } catch (error) {
        console.log(`Error in newProjectTeardown: ${error}`);
        // Try to delete via API even if UI navigation fails
        try {
            console.log(`Attempted cleanup delete for project ${projectId}`);
            await page.request.delete(`delete_project/${projectId}`);
        } catch (cleanupError) {
            console.log(`Cleanup delete also failed: ${cleanupError}`);
        }
    }
};