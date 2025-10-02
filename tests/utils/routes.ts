import type { Page } from '@playwright/test';
import type { MockProject } from './data';

//! Single source of truth for entry navigation in tests
// Update PATH to '/catalog_dev' for subpath like '/catalog_dev' for 'localhost:5170/catalog_dev' and '/' for 'localhost:5055'
const PATH = '/catalog_dev';
export async function gotoPath(page: Page) {
  await page.goto(PATH);
}

//! Update this to match the match the path of the app you are running: '/catalog_dev' for 'localhost:5170/catalog_dev' and '/' for 'localhost:5055'
export async function mockApiRoot(page: Page) {
  await page.route('**/api_root', r => r.fulfill({ json: { mdv_api_root: '/catalog_dev' } }));
  return '/catalog_dev';
}

export async function mockProjects(page: Page, items: MockProject[]) {
  await page.route('**/projects', r => r.fulfill({ json: items }));
}


