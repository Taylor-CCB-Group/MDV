# JSON Writing Stress Tests

There were intermittent issues with writing JSON files in the past (possibly exacerbated by running React with `<StrictMode>` and some fetches running twice as a result - it could be argued that it was doing a good job at highlighting the issues). This directory contains stress tests to ensure that the JSON writing process is robust and reliable.

These tests are not meant for running in CI, but as manual diagnostic tools. If we have other scenarios that we want to test for other aspects of the system, we can add them here as well.

The code was written by `@coderabbitai` [in this PR](https://github.com/Taylor-CCB-Group/MDV/pull/193) and seems fairly reasonable.

## Running the Tests

As of this writing it is hard-coded to use `"/app/mdv/json_test"` as an output directory for test files, such that it tests the writing process in the context of the volume mount. There isn't a way to change this right now, but it is a simple change to make if we need to.
To run the tests, you can use the following command:

```bash
python -m mdvtools.stress_tests.test_json_write
```

## Notes

The `"interrupted_write"` tests getting `"INVALID_JSON"` is a bit concerning, but only reflects what would happen with something crashing during a write etc so hopefully ok.

Results recorded in devcontainer during development of this test, with different strategies for synchronization (this could be parameterized, but this was tested by manually changing the code and running the tests):


### 1. No extra sync around `os.replace(temp_path, path)`:

A lot more reliable with the atomic writes, but less than 100% success rate, so we consider other strategies.

```json
{
  "interrupted_write": {
    "old": "VALID_JSON",
    "atomic": "VALID_JSON"
  },
  "concurrent_write": {
    "old": 0.2,
    "atomic": 0.95
  },
  "rapid_succession": {
    "old": 0.88,
    "atomic": 0.93
  },
  "resource_constrained": {
    "old": 0.95,
    "atomic": 1.0
  }
}
```

### 2. Using `os.system("sync")` before and after replacing the file

All of the "atomic" versions hit 100% success consistently (although the "interrupted_write" often reports "INVALID_JSON").

Potentially adds a lot of extra overhead for something that might not be reflective of real-world usage... maybe these artificially induced errors are unlikely.

```python
def save_json_atomic(path, data):
    dir_name = os.path.dirname(path)
    with tempfile.NamedTemporaryFile("w", dir=dir_name, delete=False) as tmp:
        json.dump(data, tmp, indent=2, allow_nan=False)
        tmp.flush()
        os.fsync(tmp.fileno())
        temp_name = tmp.name
    dir_name = os.path.dirname(path)
    os.system("sync")
    os.replace(temp_name, path)  # Atomic move on most OSes
    os.system("sync")

```

```json
{
  "interrupted_write": {
    "old": "VALID_JSON",
    "atomic": "INVALID_JSON"
  },
  "concurrent_write": {
    "old": 0.2,
    "atomic": 1.0
  },
  "rapid_succession": {
    "old": 0.86,
    "atomic": 1.0
  },
  "resource_constrained": {
    "old": 0.95,
    "atomic": 1.0
  }
}
```

### 3. Using `os.system("sync")` only after replacing the file

```python
def save_json_atomic(path, data):
    dir_name = os.path.dirname(path)
    with tempfile.NamedTemporaryFile("w", dir=dir_name, delete=False) as tmp:
        json.dump(data, tmp, indent=2, allow_nan=False)
        tmp.flush()
        os.fsync(tmp.fileno())
        temp_name = tmp.name
    # os.system("sync")
    os.replace(temp_name, path)  # Atomic move on most OSes
    os.system("sync")
```


```json
{
  "interrupted_write": {
    "old": "VALID_JSON",
    "atomic": "INVALID_JSON"
  },
  "concurrent_write": {
    "old": 0.2,
    "atomic": 1.0
  },
  "rapid_succession": {
    "old": 0.86,
    "atomic": 0.99
  },
  "resource_constrained": {
    "old": 0.95,
    "atomic": 1.0
  }
}
```

### 4. Using sync on the directory which should be faster

Some runs show this as having fairly significant numbers of failures, but we haven't really probably measured... maybe it is actually ok in practice, and I don't know how different the overhead is.

```python
def save_json_atomic(path, data):
    dir_name = os.path.dirname(path)
    with tempfile.NamedTemporaryFile("w", dir=dir_name, delete=False) as tmp:
        json.dump(data, tmp, indent=2, allow_nan=False)
        tmp.flush()
        os.fsync(tmp.fileno())
        temp_name = tmp.name
    dir_name = os.path.dirname(path)
    os.replace(temp_name, path)  # Atomic move on most OSes
    # this method is lower overhead than os.system("sync") 
    # but stress testing indicates it is less robust
    dir_fd = os.open(dir_name, os.O_DIRECTORY)
    os.fsync(dir_fd)
    os.close(dir_fd)
```

```json
{
  "interrupted_write": {
    "old": "VALID_JSON",
    "atomic": "INVALID_JSON"
  },
  "concurrent_write": {
    "old": 0.1,
    "atomic": 0.9
  },
  "rapid_succession": {
    "old": 0.88,
    "atomic": 0.99
  },
  "resource_constrained": {
    "old": 0.9,
    "atomic": 1.0
  }
}
```
