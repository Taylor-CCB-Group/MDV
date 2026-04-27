# Mock AnnData Module

This module provides utilities for creating realistic mock AnnData objects for testing and stress testing the MDV conversion pipeline.

## Overview

The `mock_anndata.py` module contains:

- **MockAnnDataFactory**: A factory class for creating AnnData objects with various configurations
- **Convenience functions**: Simple functions for common use cases
- **Utility functions**: For validation and summarization of AnnData objects
- **Context managers**: For suppressing warnings during testing

## Quick Start

### Basic Usage

```python
from mdvtools.tests.mock_anndata import create_minimal_anndata, create_realistic_anndata

# Create minimal AnnData for quick testing
adata_minimal = create_minimal_anndata(n_cells=20, n_genes=10)

# Create realistic AnnData with typical single-cell features
adata_realistic = create_realistic_anndata(n_cells=1000, n_genes=2000)
```

### Using the Factory

```python
from mdvtools.tests.mock_anndata import MockAnnDataFactory

# Create factory with reproducible results
factory = MockAnnDataFactory(random_seed=42)

# Create different types of AnnData
adata_minimal = factory.create_minimal(50, 25)
adata_realistic = factory.create_realistic(500, 1000)
adata_large = factory.create_large(10000, 5000)
adata_edge = factory.create_edge_cases()
```

## Available Functions

### Convenience Functions

- `create_minimal_anndata(n_cells=10, n_genes=5, add_missing=False)`: Create basic AnnData
- `create_realistic_anndata(n_cells=1000, n_genes=2000, add_missing=True)`: Create realistic single-cell data
- `create_large_anndata(n_cells=10000, n_genes=5000, add_missing=True)`: Create large datasets for stress testing
- `create_edge_case_anndata()`: Create AnnData with problematic data for edge case testing

### Factory Methods

- `create_minimal()`: Basic AnnData without extra features
- `create_realistic()`: AnnData with dimensionality reductions, layers, and unstructured data
- `create_large()`: Large AnnData with sparse matrices for stress testing
- `create_edge_cases()`: AnnData with various edge cases and problematic data
- `create_with_specific_features()`: AnnData with custom categorical features

### Utility Functions

- `get_anndata_summary(adata)`: Get comprehensive summary of AnnData properties
- `validate_anndata(adata)`: Validate AnnData structure and integrity
- `suppress_anndata_warnings()`: Context manager to suppress expected warnings

## Features

### Realistic Single-Cell Data

The mock data includes typical single-cell RNA-seq features:

**Cell Metadata (obs):**
- Cell types (T-cell, B-cell, NK-cell, etc.)
- Conditions (Control, Treatment, Disease)
- Quality metrics (total_counts, n_genes_by_counts, pct_counts_mt)
- Boolean flags (is_high_quality, is_doublet)
- Patient IDs and batch information

**Gene Metadata (var):**
- Gene types (protein_coding, lncRNA, miRNA, etc.)
- Chromosome information
- Expression statistics
- Quality flags (highly_variable, mt, ribosomal)

**Expression Data:**
- Realistic count distributions (negative binomial)
- Sparse matrices for large datasets
- Multiple layers (counts, log1p, scaled)

### Dimensionality Reductions

When creating realistic or large datasets, the following are included:

**Cell Reductions:**
- PCA (X_pca)
- UMAP (X_umap)
- t-SNE (X_tsne)

**Gene Reductions:**
- PCA (PCs in varm)

### Layers

Realistic datasets include multiple expression layers:
- `counts`: Raw count data
- `log1p`: Log-normalized data
- `scaled`: Z-score normalized data

### Unstructured Data

Realistic datasets include typical scanpy analysis results:
- Neighbors graph information
- Clustering results (Leiden)
- Differential expression results

## Edge Cases

The edge case generator creates AnnData with problematic data:

**Problematic Values:**
- Empty strings
- Very long strings
- Special characters and Unicode
- Mixed data types
- All NaN columns
- Infinity values
- Zero variance columns

**Problematic Categorical Data:**
- Empty categories
- Single category columns
- Many unique categories
- Mixed data types in categorical columns

## Manual integration testing (unique column)

To create an MDV project that includes a unique-typed column so you can test editing from the frontend and verify round-trip / `set_column_with_raw_data` behaviour:

```bash
# From repo root, with mdv environment active
python -m mdvtools.tests.example_mock_usage --create-unique-project ./test_unique_project
```

Then serve the project and open it in the browser; the `cells` datasource will have a `cell_id` unique column you can edit and save. Alternatively use `generate_test_data` with `--with-unique-column`:

```bash
python -m mdvtools.tests.generate_test_data ~/mdv/test_unique --mock --with-unique-column
```

## Testing Examples

### Basic Conversion Testing

```python
from mdvtools.tests.mock_anndata import create_minimal_anndata, suppress_anndata_warnings
from mdvtools.conversions import convert_scanpy_to_mdv

# Create test data
adata = create_minimal_anndata(100, 50)

# Test conversion
with temp_mdv_project() as test_dir:
    with suppress_anndata_warnings():
        mdv = convert_scanpy_to_mdv(test_dir, adata, delete_existing=True)
    
    # Verify conversion
    assert "cells" in mdv.get_datasource_names()
    assert "genes" in mdv.get_datasource_names()
```

### Stress Testing

```python
from mdvtools.tests.mock_anndata import MockAnnDataFactory

factory = MockAnnDataFactory(random_seed=42)

# Test with large dataset
adata_large = factory.create_large(10000, 5000)

with temp_mdv_project() as test_dir:
    with suppress_anndata_warnings():
        mdv = convert_scanpy_to_mdv(test_dir, adata_large, delete_existing=True)
    
    # Verify large dataset conversion
    cells_metadata = mdv.get_datasource_metadata("cells")
    assert len(cells_metadata['columns']) > 0
```

### Edge Case Testing

```python
from mdvtools.tests.mock_anndata import create_edge_case_anndata

# Test with problematic data
adata_edge = create_edge_case_anndata()

with temp_mdv_project() as test_dir:
    with suppress_anndata_warnings():
        mdv = convert_scanpy_to_mdv(test_dir, adata_edge, delete_existing=True)
    
    # Verify edge cases are handled gracefully
    assert isinstance(mdv, MDVProject)
```

### Custom Data Testing

```python
from mdvtools.tests.mock_anndata import MockAnnDataFactory

factory = MockAnnDataFactory(random_seed=42)

# Create data with specific features
custom_cell_types = ['Neuron', 'Astrocyte', 'Oligodendrocyte']
custom_conditions = ['Healthy', 'Alzheimer', 'Parkinson']

adata_custom = factory.create_with_specific_features(
    cell_types=custom_cell_types,
    conditions=custom_conditions,
    n_cells=200,
    n_genes=500
)

# Verify custom categories
assert set(adata_custom.obs['cell_type'].cat.categories) == set(custom_cell_types)
```

## Performance Testing

The module is designed to support performance and stress testing:

### Memory Usage Testing

```python
import psutil
import gc
from mdvtools.tests.mock_anndata import MockAnnDataFactory

factory = MockAnnDataFactory(random_seed=42)

# Monitor memory usage
gc.collect()
initial_memory = psutil.Process().memory_info().rss / 1024 / 1024

adata = factory.create_large(5000, 2000)
with temp_mdv_project() as test_dir:
    mdv = convert_scanpy_to_mdv(test_dir, adata, delete_existing=True)

gc.collect()
final_memory = psutil.Process().memory_info().rss / 1024 / 1024
memory_increase = final_memory - initial_memory

print(f"Memory increase: {memory_increase:.1f}MB")
```

### Speed Benchmarking

```python
import time
from mdvtools.tests.mock_anndata import MockAnnDataFactory

factory = MockAnnDataFactory(random_seed=42)

sizes = [(100, 50), (500, 200), (1000, 500), (2000, 1000)]

for n_cells, n_genes in sizes:
    start_time = time.time()
    adata = factory.create_realistic(n_cells, n_genes)
    
    with temp_mdv_project() as test_dir:
        conversion_start = time.time()
        mdv = convert_scanpy_to_mdv(test_dir, adata, delete_existing=True)
        conversion_time = time.time() - conversion_start
    
    cells_per_second = n_cells / conversion_time
    print(f"{n_cells}x{n_genes}: {cells_per_second:.0f} cells/s")
```

## Integration with Existing Tests

The mock module is designed to be backward compatible with existing tests. The original `create_minimal_anndata` function is preserved as a convenience function, so existing tests will continue to work without modification.

## Best Practices

1. **Use random seeds**: Set `random_seed` in the factory for reproducible tests
2. **Clean up resources**: Use context managers for temporary projects
3. **Suppress warnings**: Use `suppress_anndata_warnings()` for expected warnings
4. **Validate data**: Use `validate_anndata()` to ensure data integrity
5. **Monitor performance**: Track memory usage and conversion times for large datasets
6. **Test edge cases**: Use `create_edge_case_anndata()` to test robustness

## File Structure

```
python/mdvtools/tests/
├── mock_anndata.py              # Main mock module
├── test_categorical_data.py     # Updated tests using mock module
├── test_stress_performance.py   # Stress and performance tests
├── example_mock_usage.py        # Usage examples
└── README_mock_anndata.md       # This documentation
```

## Dependencies

The mock module requires:
- `numpy`
- `pandas`
- `scanpy`
- `scipy` (for sparse matrices)
- `psutil` (for performance testing)

## Contributing

When adding new features to the mock module:

1. Add comprehensive docstrings
2. Include type hints
3. Add tests for new functionality
4. Update this documentation
5. Ensure backward compatibility
6. Add examples for new features 