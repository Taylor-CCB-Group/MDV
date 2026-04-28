# Playwright Tests

These tests should be launched from the repository root so Playwright uses
`playwright.config.ts`.

## Quick Start

```bash
# 1. Install dependencies (if not already done)
npm install

# 2. Install Playwright browsers (one-time setup)
npx playwright install --with-deps

# 3. Start devcontainer app (if not already running)
docker compose -f docker-secrets.yml up -d

# 4. Run tests
npm run playwright-test
```

That's it! Tests run from your host machine and connect to the devcontainer app via `localhost:5055`.

For local Vite development without the devcontainer, start the app and point
Playwright at that server:

```bash
npm run dev -- --host 127.0.0.1 --port 5173
TEST_BASE_URL=http://127.0.0.1:5173 npm run playwright-test -- tests_playwright/catalog/catalog_view.spec.ts --project=chromium --reporter=list
```

The npm scripts call the local Playwright dependency directly. This avoids
`npx playwright` reaching out to the npm registry in restricted agent
environments.

## Common Workflows

### Run Tests with UI
```bash
npm run playwright-test-ui
```

### Run Specific Tests
```bash
npm run playwright-test -- tests_playwright/catalog/
npm run playwright-test -- tests_playwright/project/
```

### Run Tests in Debug Mode
```bash
npm run playwright-test -- --debug
```

## Prerequisites

- **Node.js** installed on your host machine (frontend developers already have this)
- **Devcontainer app** running on `localhost:5055`

## Folder Structure

- `catalog/` - Tests for catalog view (mocked backend)
- `project/` - Tests for project view (real backend)
- `utils/` - Test utility functions
- `test-data/` - Test data files
- `playwright.config.ts` - Playwright configuration

## Generating Test Data

Test MDV projects can be generated programmatically:

```bash
# Create project from mock data
# (requires mdv python environment to be activated - automatically the case inside container)
python -m mdvtools.tests.generate_test_data ~/mdv/test_mock --mock

# Create project from scanpy dataset
python -m mdvtools.tests.generate_test_data ~/mdv/test_pbmc3k --scanpy pbmc3k_processed

# Create large mock project
python -m mdvtools.tests.generate_test_data ~/mdv/test_large --mock --n-cells 1000000 --n-genes 5000
```

For programmatic use, import from `mdvtools.tests.test_project_factory`.

## How It Works

- **Local Execution**: Tests run on your host machine (not in a container)
- **Connect to App**: Tests connect to the devcontainer app via `localhost:5055` (exposed port)
- **Simple Setup**: Just install Playwright locally - no container configuration needed
- **UI Mode**: Works natively on your host machine

## Configuration

The `playwright.config.ts` file is located at the repository root. Key settings:

- **baseURL**: Defaults to `http://localhost:5055` (devcontainer exposed port)
- **Override**: Set `TEST_BASE_URL` environment variable to use a different URL

## Troubleshooting

### Tests can't connect to app
```bash
# Check if app service is running
docker compose -f docker-secrets.yml ps mdv_app

# Check app logs
docker compose -f docker-secrets.yml logs mdv_app

# Verify app is accessible
curl http://localhost:5055
```

### Playwright not found
```bash
# Make sure you've installed dependencies
npm install

# Install browsers
npx playwright install --with-deps
```

### Tests are failing unexpectedly
```bash
# Run with more verbose output
npm run playwright-test -- --reporter=list

# Run in headed mode to see what's happening
npm run playwright-test -- --headed
```

### Need to reset database
```bash
# Stop and remove volumes (WARNING: deletes all data)
docker compose -f docker-secrets.yml down -v
docker compose -f docker-secrets.yml up -d
```

## CI/CD

Tests run in CI using Docker services (different from local development). See `.github/workflows/playwright.yml` for the workflow configuration.

## Current Test Coverage

- **Catalog tests** (`catalog/`): Test catalog view operations with mocked backend APIs
- **Project tests** (`project/`): Test project view with real backend APIs
  - Chart creation for all chart types
  - View management operations

## Future Improvements

- Run tests in parallel (currently some run serially)
- Expand test coverage for project view
- Add more mock project generation utilities
