# UMTS (Universal Modality Translator System) - Phase 1

## Overview

Phase 1 of UMTS adds **lineage tracking** to MDVTools, providing full provenance and reproducibility for data conversions.

## What's Included

### Core Components

1. **LineageTracker** (`lineage.py`)
   - Records source files with SHA256 hashes
   - Captures conversion parameters
   - Tracks environment (Python version, packages)
   - Supports custom annotations
   - Saves/loads from JSON

2. **Utility Functions** (`utils.py`)
   - `compute_file_hash()` - SHA256 file hashing
   - `get_package_versions()` - Package version detection
   - `get_timestamp()` - ISO8601 timestamps
   - `get_file_metadata()` - File size and modification time
   - `get_environment_info()` - Full environment capture

### Integration

3. **Enhanced convert_scanpy_to_mdv()** (`conversions.py`)
   - New parameters: `track_lineage=True`, `source_file=None`
   - Automatically creates `lineage.json` in MDV projects
   - Backward compatible (existing code works unchanged)

4. **MDVProject.get_lineage()** (`mdvproject.py`)
   - Query lineage information from any MDV project
   - Returns None if no lineage available

### Testing

5. **Comprehensive Test Suite**
   - 19 unit tests for LineageTracker and utilities
   - 8 integration tests for convert_scanpy_to_mdv
   - All tests passing ✓

### Documentation

6. **User Documentation**
   - `docs/UMTS_LINEAGE.md` - Complete user guide
   - `examples/umts_lineage_demo.ipynb` - Interactive tutorial
   - Inline docstrings for all functions

## Quick Start

```python
import scanpy as sc
from mdvtools import convert_scanpy_to_mdv

# Load data
adata = sc.read_h5ad("data.h5ad")

# Convert with lineage tracking (enabled by default)
mdv = convert_scanpy_to_mdv(
    folder="./my_project",
    scanpy_object=adata,
    track_lineage=True,
    source_file="data.h5ad"
)

# Query lineage
lineage = mdv.get_lineage()
print(f"Created: {lineage['created_timestamp']}")
print(f"SHA256: {lineage['source_files'][0]['sha256']}")
```

## Files Created

```
python/mdvtools/umts/
├── __init__.py              # Module exports
├── lineage.py               # LineageTracker class
├── utils.py                 # Utility functions
└── README.md                # This file

python/mdvtools/tests/test_umts/
├── __init__.py
├── test_lineage.py          # Unit tests
└── test_integration.py      # Integration tests

docs/
└── UMTS_LINEAGE.md          # User documentation

examples/
└── umts_lineage_demo.ipynb  # Tutorial notebook
```

## Test Results

```
test_lineage.py ................ 19 passed in 0.33s ✓
test_integration.py ............  8 passed in 1.21s ✓
                                 27 tests total
```

## Benefits

- **Reproducibility**: Know exactly what went into each project
- **Auditability**: Full provenance for regulatory compliance
- **Debugging**: Track parameter and version changes
- **Documentation**: Self-documenting projects

## Future Phases (Not Implemented Yet)

- Phase 2: Global ID system for entity tracking
- Phase 3: Multi-modal integration with lineage
- Phase 4: Schema governance
- Phase 5: AI/MCP adaptive merging

## Version

UMTS Version: 0.1.0  
Phase: 1 (Lineage Tracking)  
Status: Complete ✓

