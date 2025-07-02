#!/usr/bin/env python3
"""
Test suite for categorical data handling in AnnData to MDV conversion.

This test verifies that categorical columns are properly converted and preserved.

`if data.dtype == "category"` was somewhat rashly replaced with `isinstance(data, pandas.CategoricalDtype)`
which always returned False because `data` is a pandas Series, not a dtype object.

`is_categorical_dtype(data)` properly detects categorical Series - but is deprecated for some reason
and the documentation suggested using `isinstance(dtype, pandas.CategoricalDtype)` instead.
https://pandas.pydata.org/docs/reference/api/pandas.api.types.is_categorical_dtype.html

The basic problem that caused the regression was checking `data` rather than `data.dtype`.

The motivation for the change was a type error which was caused by a pre-existing problem with the type annotation allowing for `DataFrame`.
The type of `data` passed to `add_column_to_group` should indeed always be a `Series` - but internally `add_annotations` didn't have enough
information to infer this. We now assert that internally.

Note that apparently `data = data.fillna("ND")` could also go wrong if the `ND` category is not added to the categories first.
So we also include tests to justify the assertion that this could have caused hypothetical problems before.
"""

import os
import tempfile
import shutil
import warnings
import pandas as pd
import numpy as np
import scanpy as sc
import pytest
from contextlib import contextmanager

from mdvtools.conversions import convert_scanpy_to_mdv
from mdvtools.mdvproject import MDVProject
from .mock_anndata import (
    create_minimal_anndata, 
    create_realistic_anndata,
    create_large_anndata,
    create_edge_case_anndata,
    MockAnnDataFactory,
    suppress_anndata_warnings,
    get_anndata_summary,
    validate_anndata
)


@contextmanager
def temp_mdv_project():
    """Context manager for temporary MDV project creation and cleanup."""
    test_dir = tempfile.mkdtemp()
    try:
        yield test_dir
    finally:
        if os.path.exists(test_dir):
            shutil.rmtree(test_dir)


# The create_minimal_anndata function and suppress_anndata_warnings context manager
# are now imported from the mock_anndata module


def assert_categorical_column_metadata(mdv, datasource, column_name, expect_nd=False):
    """Assert that a categorical column has correct metadata."""
    metadata = mdv.get_datasource_metadata(datasource)
    columns = {col['field']: col for col in metadata['columns']}
    
    assert column_name in columns, f"Column {column_name} not found in {datasource}"
    col_meta = columns[column_name]
    assert col_meta['datatype'] in ['text', 'text16'], f"Column {column_name} should be text type"
    assert 'values' in col_meta and len(col_meta['values']) > 0, f"Column {column_name} should have values"
    
    if expect_nd:
        assert "ND" in col_meta['values'], f"Column {column_name} should include 'ND' for missing values"


def assert_data_retrieval(mdv, datasource, column_name, expected_length, expect_nd=False):
    """Assert that data can be retrieved and has expected properties."""
    data = mdv.get_column(datasource, column_name)
    assert len(data) == expected_length, f"Data length mismatch for {column_name}"
    assert all(isinstance(val, str) for val in data), f"All values in {column_name} should be strings"
    
    if expect_nd:
        assert "ND" in data, f"'ND' should be present in {column_name} data"


def test_categorical_data_conversion():
    """Test basic categorical data conversion."""
    adata = create_minimal_anndata()
    
    with temp_mdv_project() as test_dir:
        with suppress_anndata_warnings():
            mdv = convert_scanpy_to_mdv(test_dir, adata, delete_existing=True)
        
        # Basic structure
        assert isinstance(mdv, MDVProject)
        assert "cells" in mdv.get_datasource_names()
        assert "genes" in mdv.get_datasource_names()
        
        # Test categorical columns
        for col in ['cell_type', 'condition']:
            assert_categorical_column_metadata(mdv, "cells", col)
            assert_data_retrieval(mdv, "cells", col, 10)
        
        for col in ['gene_type', 'chromosome']:
            assert_categorical_column_metadata(mdv, "genes", col)
            assert_data_retrieval(mdv, "genes", col, 5)
        
        # Test numeric columns
        metadata = mdv.get_datasource_metadata("cells")
        quality_col = next(col for col in metadata['columns'] if col['field'] == 'quality_score')
        assert quality_col['datatype'] in ['double', 'integer']


def test_category_detection_edge_cases():
    """Test the original problem: DataFrame vs Series dtype access.
    Note that this wasn't really a problem - the annotation should have been a Series.
    """
    # This demonstrates the original bug: DataFrame has no dtype attribute
    df = pd.DataFrame({'a': pd.Series(['A', 'B', 'C'], dtype='category')})
    with pytest.raises(AttributeError):
        _ = df.dtype == "category"
    
    # Series works correctly
    s = pd.Series(['A', 'B', 'C'], dtype='category')
    assert s.dtype == "category"
    
    # Current implementation should work
    assert hasattr(s, 'cat') and hasattr(s.cat, 'categories')


def test_add_annotations_series_type_assertion():
    """Test that add_annotations properly asserts Series type for individual columns."""
    # Create a minimal MDV project
    with temp_mdv_project() as test_dir:
        mdv = MDVProject(test_dir, delete_existing=True)
        
        # Create a simple datasource first
        df = pd.DataFrame({
            'id': range(5),
            'value': [1, 2, 3, 4, 5]
        })
        mdv.add_datasource("test", df)
        
        # Test that passing a DataFrame with proper structure works
        # The method expects a DataFrame with an index column and annotation columns
        annotation_df = pd.DataFrame({
            'id': [0, 1, 2, 3, 4],  # This should match the 'id' column in the datasource
            'annotation': ['A', 'B', 'C', 'D', 'E']
        })
        
        # This should work without raising an AssertionError
        mdv.add_annotations("test", annotation_df)
        
        # Verify the annotation was added
        metadata = mdv.get_datasource_metadata("test")
        columns = {col['field']: col for col in metadata['columns']}
        assert 'annotation' in columns, "Annotation column should be added"
        
        # Test that the data can be retrieved
        annotation_data = mdv.get_column("test", "annotation")
        assert len(annotation_data) == 5
        assert annotation_data[0] == 'A'
        assert annotation_data[1] == 'B'


def test_add_annotations_missing_values():
    """Test that add_annotations handles missing values correctly with the Series assert."""
    # Create a minimal MDV project
    with temp_mdv_project() as test_dir:
        mdv = MDVProject(test_dir, delete_existing=True)
        
        # Create a datasource with more rows than the annotation data
        df = pd.DataFrame({
            'id': range(10),  # 10 rows
            'value': list(range(10))
        })
        mdv.add_datasource("test", df)
        
        # Create annotation data with only some of the IDs (missing some)
        annotation_df = pd.DataFrame({
            'id': [0, 2, 4, 6, 8],  # Only even IDs
            'annotation': ['A', 'B', 'C', 'D', 'E']
        })
        
        # This should work and use the missing_value for IDs not in the annotation data
        mdv.add_annotations("test", annotation_df, missing_value="MISSING")
        
        # Verify the annotation was added
        annotation_data = mdv.get_column("test", "annotation")
        assert len(annotation_data) == 10
        
        # Check that existing IDs get their values, missing ones get the default
        assert annotation_data[0] == 'A'  # ID 0 exists
        assert annotation_data[1] == 'MISSING'  # ID 1 missing
        assert annotation_data[2] == 'B'  # ID 2 exists
        assert annotation_data[3] == 'MISSING'  # ID 3 missing


def test_categorical_missing_values():
    """Test categorical data with missing values."""
    adata = create_minimal_anndata(add_missing=True)
    
    with temp_mdv_project() as test_dir:
        with suppress_anndata_warnings():
            mdv = convert_scanpy_to_mdv(test_dir, adata, delete_existing=True)
        
        # Test that "ND" is included in categorical columns with missing values
        assert_categorical_column_metadata(mdv, "cells", "cell_type", expect_nd=True)
        assert_categorical_column_metadata(mdv, "cells", "condition", expect_nd=True)
        assert_categorical_column_metadata(mdv, "genes", "gene_type", expect_nd=True)
        
        # Test data retrieval includes "ND" values
        assert_data_retrieval(mdv, "cells", "cell_type", 10, expect_nd=True)
        assert_data_retrieval(mdv, "cells", "condition", 10, expect_nd=True)
        assert_data_retrieval(mdv, "genes", "gene_type", 5, expect_nd=True)


def test_categorical_missing_values_edge_cases():
    """Test edge cases for categorical missing value handling."""
    # Test that fillna("ND") fails without add_categories("ND")
    s = pd.Series(['A', 'B', 'C', np.nan, 'A'], dtype='category')
    
    with pytest.raises(TypeError, match="Cannot setitem on a Categorical with a new category"):
        s.fillna("ND")
    
    # Current implementation should work
    if "ND" not in s.cat.categories:
        s = s.cat.add_categories("ND")
    s = s.fillna("ND")
    
    assert "ND" in s.cat.categories
    assert s.isna().sum() == 0
    assert (s == "ND").sum() == 1


def test_boolean_dtype_handling():
    """Test boolean data conversion to categorical."""
    # Create data with boolean columns using object dtype to allow NaN
    n_cells = 10
    obs_data = pd.DataFrame({
        'is_high_quality': pd.Series([True, False, True, False, True] * 2, dtype='object'),
        'is_doublet': pd.Series([False, True, False, False, True] * 2, dtype='object'),
        'quality_score': np.random.normal(0, 1, n_cells)
    })
    
    # Add missing values
    obs_data.loc[2, 'is_high_quality'] = np.nan
    obs_data.loc[5, 'is_doublet'] = np.nan
    
    var_data = pd.DataFrame({
        'is_mitochondrial': pd.Series([False, False, True, False, False] * 2, dtype='object'),
        'mean_expression': np.random.exponential(1, 10)
    })
    
    X = np.random.negative_binomial(5, 0.3, (n_cells, 10))
    adata = sc.AnnData(X=X, obs=obs_data, var=var_data)
    
    with temp_mdv_project() as test_dir:
        with suppress_anndata_warnings():
            mdv = convert_scanpy_to_mdv(test_dir, adata, delete_existing=True)
        
        # Test boolean columns are converted to text with proper values
        for col in ['is_high_quality', 'is_doublet']:
            assert_categorical_column_metadata(mdv, "cells", col, expect_nd=True)
            assert_data_retrieval(mdv, "cells", col, n_cells, expect_nd=True)
            
            # Verify boolean values are properly converted
            data = mdv.get_column("cells", col)
            assert "True" in data and "False" in data and "ND" in data


def test_boolean_dtype_edge_cases():
    """Test edge cases for boolean conversion."""
    # Test boolean series with missing values
    s = pd.Series([True, False, True, np.nan, False], dtype='object')
    result = s.apply(lambda x: "True" if x is True else "False" if x is False else "ND")
    
    assert result[0] == "True"
    assert result[1] == "False" 
    assert result[2] == "True"
    assert result[3] == "ND"
    assert result[4] == "False"
    
    # Test all missing values
    s2 = pd.Series([np.nan, np.nan, np.nan], dtype='object')
    result2 = s2.apply(lambda x: "True" if x is True else "False" if x is False else "ND")
    assert all(val == "ND" for val in result2)


def test_categorical_missing_values_conversion_integration():
    """Test that missing value handling works in the full conversion pipeline."""
    adata = create_minimal_anndata(add_missing=True)
    
    with temp_mdv_project() as test_dir:
        with suppress_anndata_warnings():
            mdv = convert_scanpy_to_mdv(test_dir, adata, delete_existing=True)
        
        # Test specific missing value positions
        cell_type_data = mdv.get_column("cells", "cell_type")
        condition_data = mdv.get_column("cells", "condition")
        gene_type_data = mdv.get_column("genes", "gene_type")
        
        assert cell_type_data[2] == "ND", "Missing value at index 2 should be 'ND'"
        assert condition_data[5] == "ND", "Missing value at index 5 should be 'ND'"
        assert gene_type_data[3] == "ND", "Missing value at index 3 should be 'ND'"
        
        # Test non-missing values are preserved
        assert cell_type_data[0] == "T-cell"
        assert condition_data[0] == "Control"
        assert gene_type_data[0] == "protein_coding"


def test_mock_anndata_factory():
    """Test the MockAnnDataFactory class functionality."""
    factory = MockAnnDataFactory(random_seed=42)
    
    # Test minimal creation
    adata_minimal = factory.create_minimal(20, 10)
    assert adata_minimal.n_obs == 20
    assert adata_minimal.n_vars == 10
    assert not hasattr(adata_minimal, 'obsm') or len(adata_minimal.obsm) == 0
    
    # Test realistic creation
    adata_realistic = factory.create_realistic(100, 200)
    assert adata_realistic.n_obs == 100
    assert adata_realistic.n_vars == 200
    assert 'X_pca' in adata_realistic.obsm
    assert 'X_umap' in adata_realistic.obsm
    assert 'counts' in adata_realistic.layers
    assert 'neighbors' in adata_realistic.uns
    
    # Test large creation with sparse matrix
    adata_large = factory.create_large(1000, 500)
    assert adata_large.n_obs == 1000
    assert adata_large.n_vars == 500
    assert hasattr(adata_large.X, 'toarray')  # Should be sparse
    
    # Test edge cases
    adata_edge = factory.create_edge_cases()
    assert adata_edge.n_obs == 50
    assert adata_edge.n_vars == 100
    assert 'empty_string' in adata_edge.obs.columns
    assert 'all_nan' in adata_edge.obs.columns


def test_mock_anndata_with_specific_features():
    """Test creating AnnData with specific categorical features."""
    factory = MockAnnDataFactory(random_seed=42)
    
    custom_cell_types = ['Neuron', 'Astrocyte', 'Oligodendrocyte']
    custom_conditions = ['Healthy', 'Alzheimer', 'Parkinson']
    custom_gene_types = ['synaptic', 'metabolic', 'structural']
    
    adata = factory.create_with_specific_features(
        cell_types=custom_cell_types,
        conditions=custom_conditions,
        gene_types=custom_gene_types,
        n_cells=50,
        n_genes=100
    )
    
    # Check that custom categories are used
    assert set(adata.obs['cell_type'].cat.categories) == set(custom_cell_types)
    assert set(adata.obs['condition'].cat.categories) == set(custom_conditions)
    assert set(adata.var['gene_type'].cat.categories) == set(custom_gene_types)


def test_anndata_validation_and_summary():
    """Test the validation and summary utility functions."""
    factory = MockAnnDataFactory(random_seed=42)
    
    # Test with valid AnnData
    adata = factory.create_realistic(100, 200)
    assert validate_anndata(adata)
    
    summary = get_anndata_summary(adata)
    assert summary['n_cells'] == 100
    assert summary['n_genes'] == 200
    assert 'cell_type' in summary['categorical_obs']
    assert 'gene_type' in summary['categorical_var']
    assert 'X_pca' in summary['obsm_keys']
    assert 'counts' in summary['layers_keys']
    
    # Test with edge case AnnData
    adata_edge = factory.create_edge_cases()
    assert validate_anndata(adata_edge)
    summary_edge = get_anndata_summary(adata_edge)
    assert summary_edge['has_missing_obs']
    assert summary_edge['has_missing_var']


def test_stress_testing_large_datasets():
    """Test conversion with large datasets to ensure performance."""
    factory = MockAnnDataFactory(random_seed=42)
    
    # Create a moderately large dataset for testing
    adata_large = factory.create_large(5000, 2000)
    
    with temp_mdv_project() as test_dir:
        with suppress_anndata_warnings():
            mdv = convert_scanpy_to_mdv(test_dir, adata_large, delete_existing=True)
        
        # Verify conversion worked
        assert isinstance(mdv, MDVProject)
        assert "cells" in mdv.get_datasource_names()
        assert "genes" in mdv.get_datasource_names()
        
        # Check that all data was converted
        cells_metadata = mdv.get_datasource_metadata("cells")
        genes_metadata = mdv.get_datasource_metadata("genes")
        
        assert len(cells_metadata['columns']) > 0
        assert len(genes_metadata['columns']) > 0


def test_edge_case_handling():
    """Test conversion with edge case data to ensure robustness."""
    factory = MockAnnDataFactory(random_seed=42)
    
    adata_edge = factory.create_edge_cases()
    
    with temp_mdv_project() as test_dir:
        with suppress_anndata_warnings():
            # This should handle edge cases gracefully
            mdv = convert_scanpy_to_mdv(test_dir, adata_edge, delete_existing=True)
        
        # Verify conversion worked despite edge cases
        assert isinstance(mdv, MDVProject)
        assert "cells" in mdv.get_datasource_names()
        assert "genes" in mdv.get_datasource_names()


if __name__ == "__main__":
    pytest.main([__file__])