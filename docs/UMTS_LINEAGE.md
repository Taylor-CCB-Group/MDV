# UMTS Lineage Tracking

## Overview

UMTS (Universal Modality Translator System) Phase 1 adds **lineage tracking** to MDVTools, providing full provenance and reproducibility for data conversions.

When you convert single-cell data (Seurat/AnnData) to MDV format, the system now automatically records:
- Source file paths and SHA256 hashes
- Conversion parameters used
- Python and package versions
- Timestamps
- User annotations

This information is stored in a `lineage.json` file within your MDV project, creating an audit trail for regulatory compliance, debugging, and reproducibility.

## Benefits

### 1. Reproducibility
Know exactly what went into creating each MDV project. Re-run analyses with confidence that you're using the same source data and parameters.

### 2. Auditability
Track data provenance for regulatory compliance (FDA, CLIA). Every data transformation is logged with cryptographic verification.

### 3. Debugging
When results look unexpected, quickly check if the source file changed or different parameters were used.

### 4. Documentation
MDV projects are self-documenting with embedded metadata. Share projects with collaborators who can see exactly how the data was processed.

## Quick Start

### Basic Usage

```python
import scanpy as sc
from mdvtools import convert_scanpy_to_mdv

# Load your data
adata = sc.read_h5ad("pbmc3k.h5ad")

# Convert with lineage tracking (enabled by default)
mdv = convert_scanpy_to_mdv(
    folder="./pbmc3k_project",
    scanpy_object=adata,
    max_dims=3,
    track_lineage=True,          # Optional: True by default
    source_file="pbmc3k.h5ad"    # Optional: for hash tracking
)

# Query lineage information
lineage = mdv.get_lineage()
print(f"Source file: {lineage['source_files'][0]['path']}")
print(f"SHA256: {lineage['source_files'][0]['sha256']}")
print(f"Converted: {lineage['created_timestamp']}")
print(f"Python version: {lineage['environment']['python_version']}")
```

### Reopening a Project

```python
from mdvtools import MDVProject

# Open existing project
mdv = MDVProject("./pbmc3k_project")

# Query its lineage
lineage = mdv.get_lineage()

if lineage:
    print("Project created:", lineage['created_timestamp'])
    print("UMTS version:", lineage['umts_version'])
    print("Conversion function:", lineage['conversion']['function'])
    print("\nParameters used:")
    for key, value in lineage['conversion']['parameters'].items():
        print(f"  {key}: {value}")
else:
    print("No lineage information available")
```

## Lineage JSON Structure

The `lineage.json` file contains:

```json
{
  "umts_version": "0.1.0",
  "created_timestamp": "2025-11-23T10:30:00Z",
  
  "source_files": [
    {
      "path": "/data/experiments/pbmc3k.h5ad",
      "sha256": "a1b2c3d4e5f6...",
      "size_bytes": 45234567,
      "modified_timestamp": "2025-11-20T14:22:00Z"
    }
  ],
  
  "conversion": {
    "function": "convert_scanpy_to_mdv",
    "parameters": {
      "max_dims": 3,
      "delete_existing": true,
      "label": "",
      "chunk_data": false,
      "add_layer_data": true,
      "gene_identifier_column": null,
      "n_obs": 2638,
      "n_vars": 1838
    }
  },
  
  "environment": {
    "python_version": "3.10.12",
    "platform": "Linux-6.10.14-linuxkit",
    "packages": {
      "mdvtools": "1.0.0",
      "scanpy": "1.9.3",
      "scipy": "1.11.2",
      "numpy": "1.24.3",
      "pandas": "2.0.3",
      "anndata": "0.10.0",
      "h5py": "3.9.0"
    }
  },
  
  "notes": []
}
```

## Advanced Usage

### Disabling Lineage Tracking

If you need to disable lineage tracking (e.g., for testing or performance reasons):

```python
mdv = convert_scanpy_to_mdv(
    folder="./test_project",
    scanpy_object=adata,
    track_lineage=False  # Disable lineage tracking
)
```

### Adding Custom Notes

You can add notes programmatically using the LineageTracker:

```python
from mdvtools.umts import LineageTracker

# Create a tracker
tracker = LineageTracker()
tracker.record_source("input.h5ad")
tracker.record_parameters({'max_dims': 3}, function_name='convert_scanpy_to_mdv')
tracker.record_environment()

# Add custom notes
tracker.add_note("Initial conversion from Seurat export")
tracker.add_note("Filtered for quality control: removed doublets")

# Save to project
tracker.save("./my_project")
```

### Verifying File Integrity

Check if a source file has been modified since conversion:

```python
from mdvtools.umts import compute_file_hash

mdv = MDVProject("./pbmc3k_project")
lineage = mdv.get_lineage()

# Get recorded hash
original_hash = lineage['source_files'][0]['sha256']
source_path = lineage['source_files'][0]['path']

# Compute current hash
current_hash = compute_file_hash(source_path)

if current_hash == original_hash:
    print("✓ Source file unchanged")
else:
    print("⚠ Warning: Source file has been modified!")
```

## Use Cases

### Clinical Genomics Pipeline

Track patient sample processing for regulatory compliance:

```python
# Process patient sample
mdv = convert_scanpy_to_mdv(
    folder=f"./patient_{sample_id}",
    scanpy_object=adata,
    source_file=f"raw_data/{sample_id}.h5ad",
    track_lineage=True
)

# Generate audit report
lineage = mdv.get_lineage()
print("=== Audit Trail ===")
print(f"Sample processed: {lineage['created_timestamp']}")
print(f"Input file hash: {lineage['source_files'][0]['sha256']}")
print(f"MDVTools version: {lineage['environment']['packages']['mdvtools']}")
print(f"Cells analyzed: {lineage['conversion']['parameters']['n_obs']}")
```

### Cross-Study Meta-Analysis

Track which version of data was used:

```python
studies = ['study_A', 'study_B', 'study_C']

for study in studies:
    adata = sc.read_h5ad(f"{study}/processed.h5ad")
    mdv = convert_scanpy_to_mdv(
        folder=f"./meta_analysis/{study}",
        scanpy_object=adata,
        source_file=f"{study}/processed.h5ad",
        label=f"{study}_"
    )

# Later, verify all studies
for study in studies:
    mdv = MDVProject(f"./meta_analysis/{study}")
    lineage = mdv.get_lineage()
    print(f"{study}: {lineage['source_files'][0]['sha256'][:8]}...")
```

### Debugging Unexpected Results

When results change unexpectedly:

```python
# Load two versions of the same project
mdv_old = MDVProject("./backup/pbmc_project")
mdv_new = MDVProject("./current/pbmc_project")

lineage_old = mdv_old.get_lineage()
lineage_new = mdv_new.get_lineage()

# Compare source files
if lineage_old['source_files'][0]['sha256'] != lineage_new['source_files'][0]['sha256']:
    print("⚠ Different source files!")

# Compare parameters
params_old = lineage_old['conversion']['parameters']
params_new = lineage_new['conversion']['parameters']

for key in params_old:
    if params_old[key] != params_new[key]:
        print(f"Parameter changed: {key}: {params_old[key]} → {params_new[key]}")

# Compare software versions
pkg_old = lineage_old['environment']['packages']
pkg_new = lineage_new['environment']['packages']

for pkg in pkg_old:
    if pkg_old[pkg] != pkg_new.get(pkg, 'not installed'):
        print(f"Package version changed: {pkg}: {pkg_old[pkg]} → {pkg_new[pkg]}")
```

## API Reference

### LineageTracker Class

```python
from mdvtools.umts import LineageTracker

tracker = LineageTracker()
```

**Methods:**
- `record_source(file_path, compute_hash=True)` - Record a source file
- `record_parameters(params_dict, function_name=None)` - Record conversion parameters
- `record_environment()` - Record Python environment and package versions
- `add_note(message)` - Add a timestamped note
- `to_dict()` - Export lineage as dictionary
- `save(project_dir, filename='lineage.json')` - Save to file
- `load(project_dir, filename='lineage.json')` - Load from file (static method)
- `from_dict(data)` - Create from dictionary (static method)

### Utility Functions

```python
from mdvtools.umts import (
    compute_file_hash,
    get_package_versions,
    get_timestamp,
    get_file_metadata,
    get_environment_info
)
```

- `compute_file_hash(filepath, algorithm='sha256')` - Compute file hash
- `get_package_versions(packages=None)` - Get package versions
- `get_timestamp()` - Get current ISO8601 timestamp
- `get_file_metadata(filepath)` - Get file size and modification time
- `get_environment_info()` - Get Python version, platform, and packages

## Future Phases

Phase 1 (current) provides lineage tracking for single-modality conversions.

**Upcoming phases:**
- **Phase 2:** Global ID system for tracking entities across modalities
- **Phase 3:** Multi-modal integration with cross-modal lineage
- **Phase 4:** Schema governance for metadata standardization
- **Phase 5:** AI/MCP adaptive merging for intelligent data integration

## Troubleshooting

### Lineage file not created

Check that `track_lineage=True` (it's the default):

```python
mdv = convert_scanpy_to_mdv(
    folder="./my_project",
    scanpy_object=adata,
    track_lineage=True  # Explicitly enable
)
```

### SHA256 computation slow for large files

For very large files, you can skip hash computation:

```python
from mdvtools.umts import LineageTracker

tracker = LineageTracker()
tracker.record_source("large_file.h5ad", compute_hash=False)
# Hash will be recorded as "not computed"
```

### Import error

If you get an import error, ensure UMTS is installed:

```python
try:
    from mdvtools.umts import LineageTracker
    print("UMTS available")
except ImportError:
    print("UMTS not available - update mdvtools")
```

## Contributing

UMTS is part of MDVTools. To contribute:

1. Report issues on GitHub
2. Submit pull requests with tests
3. Suggest new features for future phases

## License

Same as MDVTools (see main LICENSE file).

