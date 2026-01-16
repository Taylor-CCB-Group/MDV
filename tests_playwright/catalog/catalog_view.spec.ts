import { test, expect } from '@playwright/test';
import { mockApiRoot, mockProjects, gotoPath } from '../utils/routes';
import { project } from '../utils/data';

test.describe('Catalog view', () => {
  test.beforeEach(async ({ page }) => {
    await mockApiRoot(page);
    await mockProjects(page, [project({ name: 'Alpha' }), project({ id: 'p2', name: 'Beta' })]);
    await gotoPath(page);
  });

  test('search filters', async ({ page }) => {
    await page.getByPlaceholder('Search projects').fill('Alp');
    await expect(page.getByText('Alpha')).toBeVisible();
    await expect(page.getByText('Beta')).not.toBeVisible();
  });

  //* By default, the view is grid
  test('grid/list toggle', async ({ page }) => {
    await page.getByRole('button', { name: /list view/i }).click();
    const listRows = await page.getByTestId('project_list_row').all();
    expect(listRows).toHaveLength(2);
    await page.getByRole('button', { name: /grid view/i }).click();
    const gridCards = await page.getByTestId('project_card').all();
    expect(gridCards).toHaveLength(2);
  });

 //* By default, the theme is light
  test('theme toggle', async ({ page }) => {
    const currentTheme = await page.locator('html').getAttribute('class');
    await page.getByTestId('theme_toggle_catalog').click();
    const updatedClassName = await page.locator('html').getAttribute('class');
    expect(updatedClassName).toBe(currentTheme === 'light' ? 'dark' : 'light');
  });
});


