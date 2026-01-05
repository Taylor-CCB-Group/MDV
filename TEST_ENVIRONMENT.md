# MDV Test Environment - Quick Reference

## Overview
This project uses a separate Docker container and database for Playwright testing, isolated from the development environment.

## Container Ports
- **Dev Container**: `http://localhost:5055` (your main development environment)
- **Test Container**: `http://localhost:5056` (for running tests)

Both containers can run simultaneously without interfering with each other.

## Quick Commands

### Starting/Stopping Test Environment

```bash
# Start test environment
npm run test:container:up

# Stop test environment
npm run test:container:down

# Restart with fresh database
npm run test:container:restart

# View logs
npm run test:container:logs

# Check status
./test-setup.sh status
```

### Running Tests

```bash
# Run all tests
npm run playwright-test

# Run with UI
npm run playwright-test-ui

# Complete test cycle (start → test → stop)
npm run test:e2e

# Clean run with fresh database
npm run test:e2e:clean
```

### Helper Script (Alternative)

```bash
# Start test environment
./test-setup.sh start

# Restart with fresh database
./test-setup.sh restart

# Stop test environment
./test-setup.sh stop

# View logs
./test-setup.sh logs

# Check container status
./test-setup.sh status

# Clean everything (removes all test data)
./test-setup.sh clean
```

## Typical Workflow

### Daily Testing
```bash
# 1. Start test container (do this once)
npm run test:container:up

# 2. Run tests as many times as you want
npm run playwright-test

# 3. Stop when done for the day
npm run test:container:down
```

### Fresh Start (Clean Database)
```bash
# Clean everything and start fresh
npm run test:container:restart

# Or use the all-in-one command
npm run test:e2e:clean
```

## Understanding the Setup

### Files Created
- `docker-test.yml` - Docker Compose configuration for test environment
- `test-setup.sh` - Helper script for managing test containers
- Modified `playwright.config.ts` - Now points to test container (port 5056)
- Modified `package.json` - Added convenient npm scripts

### Docker Volumes
Test environment uses separate volumes:
- `postgres-data-test` - Test database storage
- `mdv-data-test` - Test app data storage

When you run `clean` or `restart`, these volumes are deleted, giving you a fresh database.

### Database Configuration
Test container uses different database credentials:
- Database name: `mdv_test_db`
- User: `test_admin`
- Password: `test_password`
- Host: `mdv_db_test`

## Troubleshooting

### Test container won't start
```bash
# Check if port 5056 is already in use
lsof -i :5056

# View logs to see what went wrong
npm run test:container:logs
```

### Tests are failing unexpectedly
```bash
# Restart with fresh database
npm run test:container:restart

# Wait 10 seconds for services to be ready
# Then run tests again
npm run playwright-test
```

### Need to completely clean everything
```bash
# Remove all test containers and data
npm run test:container:clean

# Remove Docker images (if needed)
docker compose -f docker-test.yml down --rmi all -v
```

### Want to access test database directly
The test PostgreSQL database is exposed on port 5433:
```bash
# Connect using psql
psql -h localhost -p 5433 -U test_admin -d mdv_test_db
# Password: test_password
```

## CI/CD Integration

For GitHub Actions or other CI systems:

```yaml
# Example GitHub Actions workflow
- name: Start test environment
  run: npm run test:container:up

- name: Wait for services
  run: sleep 15

- name: Run Playwright tests
  run: npm run playwright-test

- name: Cleanup
  if: always()
  run: npm run test:container:down
```

## Tips

1. **Keep test container running** during development - no need to stop/start between test runs
2. **Use `restart`** when you want a completely fresh database
3. **Use VS Code Playwright extension** for better debugging experience
4. **Check logs** if something isn't working: `npm run test:container:logs`
5. **Both containers can run together** - develop on 5055, test on 5056

## Need Help?

Run the helper script without arguments for usage information:
```bash
./test-setup.sh
```

Or check the detailed README:
```bash
cat tests_playwright/README.md
```

