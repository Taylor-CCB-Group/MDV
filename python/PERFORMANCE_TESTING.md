# Performance Testing Setup

This document describes the performance testing infrastructure for the MDV Python tools.

## Overview

Performance tests are resource-intensive and are separated from regular CI to conserve resources. They test large dataset conversions, memory efficiency, and conversion speed benchmarks.

## Test Organization

### Regular Tests (CI)
- **Command**: `make test`
- **Scope**: All tests except those marked with `@pytest.mark.performance`
- **CI**: Runs on every commit/PR
- **Duration**: ~20 seconds

### Performance Tests (Conditional CI)
- **Command**: `make test-performance`
- **Scope**: Only tests marked with `@pytest.mark.performance`
- **CI**: Runs only when triggered (see triggers below)
- **Duration**: 2-5 minutes (depending on dataset size)

## CI Workflow

### Regular CI (`py_ci.yml`)
- Runs on every push/PR to `main` branch
- Excludes performance tests using `-m "not performance"`
- Includes pytest and pyright checks

### Performance CI (`py_performance.yml`)
- Runs only when triggered
- Works on all branches (not just main/main)
- Runs on any file changes (not just python/**)
- Includes all performance tests
- Uploads test results as artifacts

## Triggers for Performance Tests

Performance tests can be triggered in three ways:

### 1. Commit Message
Include `[perf-test]` in your commit message:
```bash
git commit -m "feat: add new conversion feature [perf-test]"
```

### 2. Pull Request Title
Include `[perf-test]` in the PR title:
```
Add new conversion feature [perf-test]
```

### 3. Manual Workflow Dispatch
1. Go to GitHub Actions tab
2. Select "Python Performance Tests" workflow
3. Click "Run workflow"
4. Choose branch and confirm

## Test Categories

### Performance Tests (`@pytest.mark.performance`)
- Large dataset conversions (20k-100k cells)
- Memory efficiency testing
- Conversion speed benchmarks
- Memory leak detection
- Chunked normalization performance

### Regular Tests
- Unit tests
- Integration tests
- Edge case handling
- Data validation

## Local mainelopment

### Running Tests Locally

```bash
# Regular tests (fast)
make test

# Performance tests (slow, resource-intensive)
make test-performance

# All tests
poetry run pytest mdvtools/
```

### Test Markers

```bash
# Run only performance tests
poetry run pytest -m performance

# Run everything except performance tests
poetry run pytest -m "not performance"

# Run slow tests
poetry run pytest -m slow

# Run integration tests
poetry run pytest -m integration
```

## Performance Test Dataset Sizes

| Test | Cells | Genes | Purpose |
|------|-------|-------|---------|
| Large | 20,000 | 3,000 | Basic performance validation |
| Very Large | 50,000 | 5,000 | Memory efficiency testing |
| Massive | 100,000 | 8,000 | Stress testing |
| Benchmarks | Variable | Variable | Performance profiling |

## Memory Requirements

Performance tests require significant memory:
- **Large tests**: ~2-4GB RAM
- **Very Large tests**: ~4-8GB RAM  
- **Massive tests**: ~8-16GB RAM

## Best Practices

1. **Use triggers sparingly**: Only run performance tests when making changes that could affect performance
2. **Monitor CI resources**: Performance tests consume significant CI minutes
3. **Local testing**: Run performance tests locally before pushing if possible
4. **Document changes**: Include performance impact in commit messages when relevant

## Troubleshooting

### Workflow Not Triggering
If performance tests are not running despite having `[perf-test]` in your commit message:

1. **Check your branch**: The workflow now runs on all branches (not just main/main)
2. **Verify commit message**: Make sure `[perf-test]` is exactly as shown (case-sensitive)
3. **Check file changes**: The workflow runs on any file changes (not just python files)
4. **Debug locally**: Check the workflow logs for debugging information

### Tests Getting Killed
If performance tests are killed due to memory:
1. Reduce dataset sizes in test configuration
2. Run tests on a machine with more RAM
3. Use CI with higher resource limits

### CI Timeouts
If CI times out:
1. Reduce the number of benchmark iterations
2. Use smaller test datasets
3. Split tests into smaller chunks

## Configuration

Test configuration is in `pyproject.toml`:
```toml
[tool.pytest.ini_options]
markers = [
    "performance: marks tests as performance tests",
    "slow: marks tests as slow running",
    "integration: marks tests as integration tests"
]
``` 