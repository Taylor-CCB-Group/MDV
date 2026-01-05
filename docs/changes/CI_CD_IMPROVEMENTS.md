# CI/CD Pipeline Improvements

This document outlines the recent refinements to the CI/CD pipeline, aimed at improving security, clarity, and build versioning.

## 1. Consolidation of Workflows

Previously, the repository contained multiple, overlapping GitHub Actions workflows (`deploy.yml`, `daily-build-check.yml`, `deploy_with_labels.yml`), which caused confusion.

**Change:** These have been consolidated into a single, unified workflow file: `.github/workflows/build-and-deploy.yml`.

This new workflow handles all build and deployment logic and is triggered by:
- A daily schedule at midnight.
- A push to the `main` branch with `[deploy]` in the commit message.
- A pull request being merged into the `main` branch.
- Manual dispatch via the GitHub Actions UI.

## 2. Enhanced Security

The previous workflow required write permissions to the repository to commit version changes, posing a potential security risk.

**Change:** The need to commit from within the CI pipeline has been eliminated.
- The `permissions: contents: write` block has been removed from the workflow.
- The `GITHUB_TOKEN` now uses the default, read-only permissions for builds, adhering to the principle of least privilege.

## 3. Improved Build Versioning and Metadata

The old system relied on a simple, manually incremented `version.txt` file, which provided limited traceability.

**Change:** We have moved to a more robust versioning strategy using Git metadata.
- **Removed `version.txt`:** This file is no longer used.
- **Git Commit SHA as Version:** Docker images are now tagged with the Git commit SHA (e.g., `your-repo/your-image:a1b2c3d`), creating a direct, immutable link between the image and the source code it was built from.
- **Rich Build Information:** The following build-time metadata is now injected directly into the Docker image as environment variables:
    - Commit SHA (`VITE_GIT_COMMIT_HASH`)
    - Commit Date (`VITE_GIT_COMMIT_DATE`)
    - Branch Name (`VITE_GIT_BRANCH_NAME`)
    - Last Commit Message (`VITE_GIT_LAST_COMMIT_MESSAGE`)
    - Build Timestamp (`VITE_BUILD_DATE`)
    - Git Dirty Status (for local builds) (`VITE_GIT_DIRTY`)

## 4. Consistent Local and CI Environments

It was difficult to get consistent build information between local development and CI builds.

**Change:** The system is now aligned to provide the same metadata in all environments.
- **`vite.config.mts`:** Updated to gather local Git information when running the dev server.
- **`Dockerfile`:** Updated to accept build arguments for all the metadata listed above.
- **`build-and-deploy.yml`:** The CI workflow now generates this metadata and passes it to the `docker build` command.

This ensures that the `useBuildInfo()` hook in the frontend receives the same data structure, whether running locally or from a deployed container.

## 5. Local development and Building

To ensure build information is consistent during development, the following approach is used:

- **Live development:** When running the Vite dev server inside the devcontainer (`npm run dev`), the `vite.config.mts` file automatically detects the mounted `.git` directory and injects the current Git status into the application.
- **Building a Local Image:** To build a production-like Docker image locally that contains the correct Git metadata, the configuration is handled by `docker-compose`. The `docker-secrets.yml` file (used by the dev container) is configured to pass build arguments from the host environment to the Docker build process. You can build the image by running `docker compose -f docker-secrets.yml build` in your terminal, which will automatically pass the necessary environment variables. 