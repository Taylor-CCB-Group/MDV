import { test, expect } from '@playwright/test';
import { mockApiRoot, mockProjects, gotoPath } from '../utils/routes';

test('create project navigates to project view', async ({ page }) => {
  await mockApiRoot(page);
  await mockProjects(page, []);

  // mock create project api
  await page.route('**/create_project', r => r.fulfill({ json: { id: '123', name: 'New Project' } }));

  await gotoPath(page);
  await page.getByText('Create new project').click();
  await page.waitForURL('**/project/123');
  expect(page.url()).toContain('/project/123');
});