import { expect, test } from "@playwright/test";
import { mockApiRoot, mockProjects, gotoPath } from "../utils/routes";
import { project } from "../utils/data";

test.describe('share project', () => {
  test.beforeEach(async ({ page }) => {
    await mockApiRoot(page);
    await mockProjects(page, [project({ id: 'p1', name: 'Alpha' })]);
    await gotoPath(page);
  });

  test('add new user', async ({ page }) => {
    await page.route('**/projects/p1/share', async (route) => {
      if (route.request().method() !== 'GET') {
        await route.fallback();
        return;
      }
      await route.fulfill({ json: { all_users: [{ id: 1, email: 'test@test.com' }, { id: 2, email: 'test2@test.com' }], shared_users: [{ id: 1, email: 'test@test.com', permission: 'View' }] } });
    });
    await page.route('**/projects/p1/share', async (route) => {
      if (route.request().method() !== 'POST') {
        await route.fallback();
        return;
      }
      await route.fulfill({ json: { message: 'ok' } });
    });
    await page.getByTestId('project_menu_p1').click();
    await page.getByTestId('project_share_p1').click();
    await expect(page.getByText(/manage project sharing/i)).toBeVisible();
    await expect(page.getByText('test@test.com')).toBeVisible();
    await page.getByPlaceholder('Enter email to search for the user').fill('test2@test.com');
    await page.getByRole('option', { name: 'test2@test.com' }).click();
    await page.route('**/projects/p1/share', async (route) => {
      if (route.request().method() !== 'GET') {
        await route.fallback();
        return;
      }
      await route.fulfill({ json: { all_users: [{ id: 1, email: 'test@test.com' }, { id: 2, email: 'test2@test.com' }], shared_users: [{ id: 1, email: 'test@test.com', permission: 'View' }, { id: 2, email: 'test2@test.com', permission: 'View' }] } });
    });
    await page.getByRole('button', { name: /add/i }).click();
    await expect(page.getByText('test2@test.com')).toBeVisible();
  });

  // todo: add locators for accessing the user permission dropdown
  test('change user permission', async ({ page }) => {
    await page.route('**/projects/p1/share/2/edit', async (route) => {
      if (route.request().method() !== 'POST') {
        await route.fallback();
        return;
      }
      await route.fulfill({ json: { message: 'ok' } });
    });
  });

  // todo: add locators for accessing the delete button
  test('delete user', async ({ page }) => {
    await page.route('**/projects/p1/share/2/delete', async (route) => {
      if (route.request().method() !== 'POST') {
        await route.fallback();
        return;
      }
      await route.fulfill({ json: { message: 'ok' } });
    });
  });
});
