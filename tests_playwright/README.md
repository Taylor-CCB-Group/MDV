## Playwright tests

### Folder structure
- The playwright tests are located in the tests folder.
- Currently we have `/catalog` and `/project` for testing catalog and project related operations.
- Test data is located in `/test-data`
- The test utility functions are in `/utils`

### Running the tests

#### New Approach: Isolated Test Environment (Recommended)
We now have a separate Docker container and database specifically for testing. This provides:
- ✅ **Isolated environment**: Tests don't affect your dev database
- ✅ **Fresh database**: Each test run can start with a clean state
- ✅ **Parallel testing**: No conflicts between test runs
- ✅ **CI/CD ready**: Easy to integrate into automated pipelines

**Quick Start:**
```bash
# 1. Start the test environment (runs on port 5056)
npm run test:container:up

# 2. Run tests
npm run playwright-test
# OR use VS Code Playwright extension

# 3. Stop when done (optional)
npm run test:container:down

# 4. Clean everything for fresh start (removes all test data)
npm run test:container:clean
```

**Available Commands:**
- `npm run test:container:up` - Start test containers
- `npm run test:container:down` - Stop test containers
- `npm run test:container:clean` - Stop and remove all test data
- `npm run test:container:restart` - Quick restart with fresh database
- `npm run test:container:logs` - View container logs
- `npm run test:e2e` - Run complete test cycle (start → test → stop)
- `npm run test:e2e:clean` - Clean run with fresh database

**Using the helper script (alternative):**
```bash
./test-setup.sh start    # Start test environment
./test-setup.sh restart  # Restart with fresh database
./test-setup.sh stop     # Stop test environment
./test-setup.sh logs     # View logs
./test-setup.sh status   # Check status
```

**Port Configuration:**
- Dev container: `http://localhost:5055`
- Test container: `http://localhost:5056`
- Both can run simultaneously without conflicts

#### Old Approach (Still supported)
- Before running the tests, open the dev container to run the app, and open the project outside the container in another window to run the tests.
- To run the tests, we can either run the `playwright-test` script or run them using the vscode playwright extension.
- The vscode playwright extension is recommended for running tests as it is easy and provides more features.
- Run the catalog and project folder tests separately as there was some weird behaviour of tests failing when done simultaneously.

### Current progress
- The tests in `/catalog` contains the tests related to the catalog view. It consists of tests to test the operations in the catalog view by mocking the backend APIs.
- The tests in `/project` contains the tests related to the project view. Here the backend is not mocked, we use the actual backend APIs in the tests. It consists of tests to:
    - Create all the chart types and assert if they are created
    - Perform all the view related operations
- Currently, the tests are run in serial mode in `/project` because we need to have access to a new project for all tests.
- There is no CI/CD for the playwright tests to be run. The main blocker with this was that we need a dev container to be running for the playwright tests to run, so in the CI/CD tests, there should be a dev container running.

### Future plan
- The CI/CD for playwright tests should be implemented either by having a dev container running while running the tests or by a different approach.
- The tests should run in parallel rather than serial mode. This maybe fixed using project fixtures or setup and teardown.
- Add more tests for the project view.
- Have a better way of running tests by using a separate container.