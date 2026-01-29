import { test, expect } from '@playwright/test';
import { mockApiRoot, mockProjects, gotoPath } from '../utils/routes';
import path from 'node:path';

// This test currently only works with HTTP upload, we need to update the USE_SOCKETIO_UPLOAD flg in the ImportProjectDialog.tsx to false to use HTTP upload.
// todo: update the test to work with socketio upload
test('import dialog opens and success navigates', async ({ page }) => {
  await mockApiRoot(page);
  await mockProjects(page, []);
  await gotoPath(page);
  await page.waitForLoadState('networkidle');

  await page.getByText('Import an existing project').click();
  
  // Wait for dialog to open
  await page.waitForTimeout(500);
  
  // Upload zip file
  const filePath = path.join(__dirname, '..', 'test-data', 'pbmc3k-mdv.mdv.zip');
  await page.locator('input[type="file"]').setInputFiles(filePath);

  await page.route('**/import_project', r => r.fulfill({ json: { status: 'success', id: '555' } }));

  await page.getByRole('button', { name: /upload file/i }).click();

  // Increased timeout for navigation
  await page.waitForURL('**/project/**', { timeout: 60000 });
  expect(page.url()).toContain('/project/');
});


