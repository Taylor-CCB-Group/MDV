## Playwright tests

### Folder structure
- The playwright tests are located in the tests folder.
- Currently we have `/catalog` and `/project` for testing catalog and project related operations.
- Test data is located in `/test-data`
- The test utility functions are in `/utils`

### Running the tests
- Before running the tests, ```open the dev container to run the app, and open the project outside the container in another window to run the tests```. (This is required because we need the app running to access the app, which can be done via dev container at the moment)
- To run the tests, we can either run the `playwright-test` script or run them using the vscode playwright extension.
- The vscode playwright extension is recommended for running tests as it is easy and provides more features.

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