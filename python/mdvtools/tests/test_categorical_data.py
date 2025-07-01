#!/usr/bin/env python3
"""
Test suite for categorical data handling in AnnData to MDV conversion.

This test verifies that categorical columns are properly converted and preserved.

`if data.dtype == "category"` was somewhat rashly replaced with `isinstance(data, pandas.CategoricalDtype)`
which always returned False because `data` is a pandas Series, not a dtype object.

`is_categorical_dtype(data)` properly detects categorical Series - but is deprecated for some reason
and the documentation suggested using `isinstance(dtype, pandas.CategoricalDtype)` instead.
https://pandas.pydata.org/docs/reference/api/pandas.api.types.is_categorical_dtype.html

The problem was that it was checking `data` rather than `data.dtype`.

There may be scenarios in which the original `data.dtype == "category"` approach would fail,
perhaps not likely in practice, but this test aims to cover all edge cases and verify that the
current implementation is robust and any regressions are caught.

This has mostly been written in Cursor - rather more verbose than it needs to be.

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

from mdvtools.conversions import convert_scanpy_to_mdv
from mdvtools.mdvproject import MDVProject


def create_test_anndata():
    """Create a synthetic AnnData object with categorical and numeric data for testing."""
    n_cells, n_genes = 50, 25
    
    # Cell metadata with categorical and numeric columns
    obs_data = pd.DataFrame({
        'cell_type': pd.Categorical(['T-cell', 'B-cell', 'NK-cell'] * 16 + ['T-cell', 'B-cell']),  # 50 elements
        'condition': pd.Categorical(['Control', 'Treatment'] * 25),  # 50 elements
        'quality_score': np.random.normal(0, 1, n_cells),  # 50 elements
        'is_high_quality': pd.Categorical(['Yes', 'No'] * 25)  # 50 elements
    })
    
    # Gene metadata with categorical and numeric columns
    var_data = pd.DataFrame({
        'gene_type': pd.Categorical(['protein_coding', 'lncRNA', 'pseudogene'] * 8 + ['miRNA']),  # 25 elements
        'chromosome': pd.Categorical(['chr1', 'chr2', 'chrX'] * 8 + ['chrY']),  # 25 elements
        'mean_expression': np.random.exponential(1, n_genes),  # 25 elements
        'is_mitochondrial': pd.Categorical(['Yes', 'No'] * 12 + ['No'])  # 25 elements
    })
    
    # Expression matrix
    X = np.random.negative_binomial(5, 0.3, (n_cells, n_genes))
    
    return sc.AnnData(X=X, obs=obs_data, var=var_data)


def test_categorical_data_conversion():
    """Test that categorical data is properly converted from AnnData to MDV format."""
    adata = create_test_anndata()
    test_project_dir = tempfile.mkdtemp()
    
    try:
        # Convert to MDV format - suppress expected AnnData warnings
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="Transforming to str index", category=UserWarning)
            mdv = convert_scanpy_to_mdv(test_project_dir, adata, delete_existing=True)
        
        # Basic assertions
        assert isinstance(mdv, MDVProject)
        assert "cells" in mdv.get_datasource_names()
        assert "genes" in mdv.get_datasource_names()
        
        # Get datasource metadata
        cell_datasource = mdv.get_datasource_metadata("cells")
        gene_datasource = mdv.get_datasource_metadata("genes")
        
        # Verify categorical columns are properly handled
        cell_columns = {col['field']: col for col in cell_datasource['columns']}
        gene_columns = {col['field']: col for col in gene_datasource['columns']}
        
        # Test categorical columns in cells
        categorical_cell_cols = ['cell_type', 'condition', 'is_high_quality']
        for col_name in categorical_cell_cols:
            assert col_name in cell_columns
            col_meta = cell_columns[col_name]
            assert col_meta['datatype'] in ['text', 'text16']
            assert 'values' in col_meta and len(col_meta['values']) > 0
        
        # Test numeric columns in cells
        assert 'quality_score' in cell_columns
        assert cell_columns['quality_score']['datatype'] in ['double', 'integer']
        
        # Test categorical columns in genes
        categorical_gene_cols = ['gene_type', 'chromosome', 'is_mitochondrial']
        for col_name in categorical_gene_cols:
            assert col_name in gene_columns
            col_meta = gene_columns[col_name]
            assert col_meta['datatype'] in ['text', 'text16']
            assert 'values' in col_meta and len(col_meta['values']) > 0
        
        # Test numeric columns in genes
        assert 'mean_expression' in gene_columns
        assert gene_columns['mean_expression']['datatype'] in ['double', 'integer']
        
        # Test data retrieval and value preservation
        for col_name in categorical_cell_cols:
            retrieved_data = mdv.get_column("cells", col_name)
            assert len(retrieved_data) == 50
            assert all(isinstance(val, str) for val in retrieved_data)
            
            # Check original values are preserved
            original_unique = set(adata.obs[col_name].cat.categories)
            retrieved_unique = set(retrieved_data)
            assert original_unique.issubset(retrieved_unique)
        
    finally:
        if os.path.exists(test_project_dir):
            shutil.rmtree(test_project_dir)


def test_category_detection_edge_cases():
    """Test edge cases for categorical detection."""
    # Test DataFrame (should raise AttributeError with old approach)
    df = pd.DataFrame({'a': pd.Series(['A', 'B', 'C'], dtype='category')})
    with pytest.raises(AttributeError):
        _ = df.dtype == "category"
    
    # Test Series (works correctly)
    s = pd.Series(['A', 'B', 'C'], dtype='category')
    assert s.dtype == "category"
    
    # Test current implementation is robust
    assert hasattr(s, 'cat') and hasattr(s.cat, 'categories')


def test_categorical_missing_values():
    """Test that categorical data with missing values is handled correctly."""
    # Create test data with missing values
    n_cells = 20
    
    # Test case 1: Categorical data with NaN values
    obs_data = pd.DataFrame({
        'cell_type': pd.Categorical(['T-cell', 'B-cell', 'NK-cell', 'T-cell', 'B-cell'] * 4),
        'condition': pd.Categorical(['Control', 'Treatment', 'Control', 'Treatment', 'Control'] * 4),
        'quality_score': np.random.normal(0, 1, n_cells),
        'is_high_quality': pd.Categorical(['Yes', 'No', 'Yes', 'No', 'Yes'] * 4)
    })
    
    # Introduce missing values
    obs_data.loc[2, 'cell_type'] = np.nan
    obs_data.loc[5, 'condition'] = np.nan
    obs_data.loc[8, 'is_high_quality'] = np.nan
    
    # Create gene data with missing values
    var_data = pd.DataFrame({
        'gene_type': pd.Categorical(['protein_coding', 'lncRNA', 'pseudogene'] * 8 + ['miRNA']),
        'chromosome': pd.Categorical(['chr1', 'chr2', 'chrX'] * 8 + ['chrY']),
        'mean_expression': np.random.exponential(1, 25),
        'is_mitochondrial': pd.Categorical(['Yes', 'No'] * 12 + ['No'])
    })
    
    # Introduce missing values in gene data
    var_data.loc[3, 'gene_type'] = np.nan
    var_data.loc[7, 'chromosome'] = np.nan
    
    # Create expression matrix
    X = np.random.negative_binomial(5, 0.3, (n_cells, 25))
    
    # Create AnnData object
    adata = sc.AnnData(X=X, obs=obs_data, var=var_data)
    
    test_project_dir = tempfile.mkdtemp()
    
    try:
        # Convert to MDV format - suppress expected AnnData warnings
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="Transforming to str index", category=UserWarning)
            mdv = convert_scanpy_to_mdv(test_project_dir, adata, delete_existing=True)
        
        # Verify the conversion succeeded
        assert isinstance(mdv, MDVProject)
        assert "cells" in mdv.get_datasource_names()
        assert "genes" in mdv.get_datasource_names()
        
        # Get datasource metadata
        cell_datasource = mdv.get_datasource_metadata("cells")
        gene_datasource = mdv.get_datasource_metadata("genes")
        
        # Check that categorical columns with missing values are handled
        cell_columns = {col['field']: col for col in cell_datasource['columns']}
        gene_columns = {col['field']: col for col in gene_datasource['columns']}
        
        # Test that "ND" is in the values list for columns with missing data
        categorical_cell_cols_with_missing = ['cell_type', 'condition', 'is_high_quality']
        for col_name in categorical_cell_cols_with_missing:
            assert col_name in cell_columns
            col_meta = cell_columns[col_name]
            assert col_meta['datatype'] in ['text', 'text16']
            assert 'values' in col_meta and len(col_meta['values']) > 0
            # Check that "ND" is included in the values (for missing data handling)
            assert "ND" in col_meta['values'], f"'ND' should be in values for {col_name}"
        
        # Test that "ND" is in the values list for gene columns with missing data
        categorical_gene_cols_with_missing = ['gene_type', 'chromosome']
        for col_name in categorical_gene_cols_with_missing:
            assert col_name in gene_columns
            col_meta = gene_columns[col_name]
            assert col_meta['datatype'] in ['text', 'text16']
            assert 'values' in col_meta and len(col_meta['values']) > 0
            # Check that "ND" is included in the values (for missing data handling)
            assert "ND" in col_meta['values'], f"'ND' should be in values for {col_name}"
        
        # Test data retrieval and verify missing values are properly handled
        for col_name in categorical_cell_cols_with_missing:
            retrieved_data = mdv.get_column("cells", col_name)
            assert len(retrieved_data) == n_cells
            assert all(isinstance(val, str) for val in retrieved_data)
            
            # Check that "ND" values are present in the retrieved data
            assert "ND" in retrieved_data, f"'ND' should be present in retrieved data for {col_name}"
            
            # Check that original values are preserved
            original_unique = set(adata.obs[col_name].cat.categories)
            retrieved_unique = set(retrieved_data)
            assert original_unique.issubset(retrieved_unique)
        
    finally:
        if os.path.exists(test_project_dir):
            shutil.rmtree(test_project_dir)


def test_categorical_missing_values_edge_cases():
    """Test edge cases for categorical missing value handling."""
    # Test case 1: Verify that fillna("ND") fails without add_categories("ND")
    s = pd.Series(['A', 'B', 'C', np.nan, 'A'], dtype='category')
    
    # This should raise a TypeError because "ND" is not in the categories
    with pytest.raises(TypeError, match="Cannot setitem on a Categorical with a new category"):
        s.fillna("ND")
    
    # This should work (simulating the current implementation)
    if "ND" not in s.cat.categories:
        s = s.cat.add_categories("ND")
    s = s.fillna("ND")
    
    # Verify the result
    assert "ND" in s.cat.categories
    assert s.isna().sum() == 0  # No more NaN values
    assert (s == "ND").sum() == 1  # One "ND" value
    
    # Test case 2: All missing values
    s2 = pd.Series([np.nan, np.nan, np.nan], dtype='category')
    
    # This should also fail without add_categories
    with pytest.raises(TypeError, match="Cannot setitem on a Categorical with a new category"):
        s2.fillna("ND")
    
    # This should work
    if "ND" not in s2.cat.categories:
        s2 = s2.cat.add_categories("ND")
    s2 = s2.fillna("ND")
    
    # Verify the result
    assert "ND" in s2.cat.categories
    assert s2.isna().sum() == 0
    assert (s2 == "ND").sum() == 3  # All values are "ND"
    
    # Test case 3: Empty categorical with missing values
    s3 = pd.Series([], dtype='category')
    s3 = s3.cat.add_categories("ND")
    s3 = s3.fillna("ND")
    
    # Verify the result
    assert "ND" in s3.cat.categories
    assert len(s3) == 0


def test_boolean_dtype_handling():
    """Test that boolean dtype data is handled correctly without category issues."""
    # Create test data with boolean columns
    n_cells = 20
    
    # Use object dtype to allow NaN values in boolean columns
    obs_data = pd.DataFrame({
        'is_high_quality': pd.Series([True, False, True, False, True] * 4, dtype='object'),
        'is_doublet': pd.Series([False, True, False, False, True] * 4, dtype='object'),
        'has_mitochondrial_genes': pd.Series([True, True, False, True, False] * 4, dtype='object'),
        'quality_score': np.random.normal(0, 1, n_cells)  # numeric column
    })
    
    # Introduce missing values in boolean columns
    obs_data.loc[2, 'is_high_quality'] = np.nan
    obs_data.loc[5, 'is_doublet'] = np.nan
    obs_data.loc[8, 'has_mitochondrial_genes'] = np.nan
    
    # Create gene data with boolean columns
    var_data = pd.DataFrame({
        'is_mitochondrial': pd.Series([False, False, True, False, False] * 5, dtype='object'),
        'is_ribosomal': pd.Series([True, False, False, True, False] * 5, dtype='object'),
        'is_highly_variable': pd.Series([True, False, True, False, True] * 5, dtype='object'),
        'mean_expression': np.random.exponential(1, 25)  # numeric column
    })
    
    # Introduce missing values in gene boolean columns
    var_data.loc[3, 'is_mitochondrial'] = np.nan
    var_data.loc[7, 'is_ribosomal'] = np.nan
    var_data.loc[10, 'is_highly_variable'] = np.nan  # Add missing value to this column too
    
    # Create expression matrix
    X = np.random.negative_binomial(5, 0.3, (n_cells, 25))
    
    # Create AnnData object
    adata = sc.AnnData(X=X, obs=obs_data, var=var_data)
    
    test_project_dir = tempfile.mkdtemp()
    
    try:
        # Convert to MDV format - suppress expected AnnData warnings
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="Transforming to str index", category=UserWarning)
            mdv = convert_scanpy_to_mdv(test_project_dir, adata, delete_existing=True)
        
        # Verify the conversion succeeded
        assert isinstance(mdv, MDVProject)
        assert "cells" in mdv.get_datasource_names()
        assert "genes" in mdv.get_datasource_names()
        
        # Get datasource metadata
        cell_datasource = mdv.get_datasource_metadata("cells")
        gene_datasource = mdv.get_datasource_metadata("genes")
        
        # Check that boolean columns are handled correctly
        cell_columns = {col['field']: col for col in cell_datasource['columns']}
        gene_columns = {col['field']: col for col in gene_datasource['columns']}
        
        # Test boolean columns in cells
        boolean_cell_cols = ['is_high_quality', 'is_doublet', 'has_mitochondrial_genes']
        for col_name in boolean_cell_cols:
            assert col_name in cell_columns
            col_meta = cell_columns[col_name]
            assert col_meta['datatype'] in ['text', 'text16']
            assert 'values' in col_meta and len(col_meta['values']) > 0
            # Check that boolean values are properly converted to strings
            assert "True" in col_meta['values'], f"'True' should be in values for {col_name}"
            assert "False" in col_meta['values'], f"'False' should be in values for {col_name}"
            # Check that "ND" is included for missing values
            assert "ND" in col_meta['values'], f"'ND' should be in values for {col_name}"
        
        # Test boolean columns in genes
        boolean_gene_cols = ['is_mitochondrial', 'is_ribosomal', 'is_highly_variable']
        for col_name in boolean_gene_cols:
            assert col_name in gene_columns
            col_meta = gene_columns[col_name]
            assert col_meta['datatype'] in ['text', 'text16']
            assert 'values' in col_meta and len(col_meta['values']) > 0
            # Check that boolean values are properly converted to strings
            assert "True" in col_meta['values'], f"'True' should be in values for {col_name}"
            assert "False" in col_meta['values'], f"'False' should be in values for {col_name}"
            # Check that "ND" is included for missing values
            assert "ND" in col_meta['values'], f"'ND' should be in values for {col_name}"
        
        # Test data retrieval and verify boolean values are properly handled
        for col_name in boolean_cell_cols:
            retrieved_data = mdv.get_column("cells", col_name)
            assert len(retrieved_data) == n_cells
            assert all(isinstance(val, str) for val in retrieved_data)
            
            # Check that boolean values are properly converted
            assert "True" in retrieved_data, f"'True' should be present in retrieved data for {col_name}"
            assert "False" in retrieved_data, f"'False' should be present in retrieved data for {col_name}"
            assert "ND" in retrieved_data, f"'ND' should be present in retrieved data for {col_name}"
            
            # Verify the conversion logic works correctly
            true_count = sum(1 for val in retrieved_data if val == "True")
            false_count = sum(1 for val in retrieved_data if val == "False")
            nd_count = sum(1 for val in retrieved_data if val == "ND")
            
            # Should have some of each type
            assert true_count > 0, f"Should have some 'True' values in {col_name}"
            assert false_count > 0, f"Should have some 'False' values in {col_name}"
            assert nd_count > 0, f"Should have some 'ND' values in {col_name}"
        
    finally:
        if os.path.exists(test_project_dir):
            shutil.rmtree(test_project_dir)


def test_boolean_dtype_edge_cases():
    """Test edge cases for boolean dtype handling."""
    # Test case 1: Boolean series with missing values (using object dtype)
    s = pd.Series([True, False, True, np.nan, False], dtype='object')
    
    # The current implementation should handle this correctly
    # by converting to strings and handling NaN as "ND"
    result = s.apply(lambda x: "True" if x is True else "False" if x is False else "ND")
    
    # Verify the result
    assert result[0] == "True"
    assert result[1] == "False"
    assert result[2] == "True"
    assert result[3] == "ND"
    assert result[4] == "False"
    
    # Test case 2: All missing boolean values
    s2 = pd.Series([np.nan, np.nan, np.nan], dtype='object')
    result2 = s2.apply(lambda x: "True" if x is True else "False" if x is False else "ND")
    
    # Verify the result
    assert all(val == "ND" for val in result2)
    
    # Test case 3: Empty boolean series
    s3 = pd.Series([], dtype='object')
    result3 = s3.apply(lambda x: "True" if x is True else "False" if x is False else "ND")
    
    # Verify the result
    assert len(result3) == 0
    
    # Test case 4: Pure boolean series (no missing values)
    s4 = pd.Series([True, False, True, False], dtype='bool')
    result4 = s4.apply(lambda x: "True" if x is True else "False" if x is False else "ND")
    
    # Verify the result
    assert result4[0] == "True"
    assert result4[1] == "False"
    assert result4[2] == "True"
    assert result4[3] == "False"


def test_categorical_missing_values_conversion_integration():
    """Test that the missing value handling works in the full conversion pipeline."""
    # Create a minimal test case that would fail without the add_categories logic
    n_cells = 10
    
    # Create data with missing values in categorical columns
    obs_data = pd.DataFrame({
        'cell_type': pd.Categorical(['T-cell', 'B-cell', 'NK-cell', 'T-cell', 'B-cell'] * 2),
        'condition': pd.Categorical(['Control', 'Treatment'] * 5)
    })
    
    # Introduce missing values
    obs_data.loc[2, 'cell_type'] = np.nan
    obs_data.loc[5, 'condition'] = np.nan
    
    # Create minimal gene data
    var_data = pd.DataFrame({
        'gene_type': pd.Categorical(['protein_coding', 'lncRNA', 'pseudogene', 'miRNA', 'protein_coding'] * 5),
        'chromosome': pd.Categorical(['chr1', 'chr2', 'chrX', 'chrY', 'chr1'] * 5)
    })
    
    # Introduce missing values in gene data
    var_data.loc[3, 'gene_type'] = np.nan
    var_data.loc[7, 'chromosome'] = np.nan
    
    # Create minimal expression matrix
    X = np.random.negative_binomial(5, 0.3, (n_cells, 25))
    
    # Create AnnData object
    adata = sc.AnnData(X=X, obs=obs_data, var=var_data)
    
    test_project_dir = tempfile.mkdtemp()
    
    try:
        # This conversion should succeed without errors - suppress expected AnnData warnings
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="Transforming to str index", category=UserWarning)
            mdv = convert_scanpy_to_mdv(test_project_dir, adata, delete_existing=True)
        
        # Verify conversion succeeded
        assert isinstance(mdv, MDVProject)
        
        # Test that we can retrieve data with missing values
        cell_type_data = mdv.get_column("cells", "cell_type")
        condition_data = mdv.get_column("cells", "condition")
        
        # Check that "ND" values are present where we had NaN
        assert cell_type_data[2] == "ND", "Missing value at index 2 should be 'ND'"
        assert condition_data[5] == "ND", "Missing value at index 5 should be 'ND'"
        
        # Check that non-missing values are preserved
        assert cell_type_data[0] == "T-cell"
        assert cell_type_data[1] == "B-cell"
        assert condition_data[0] == "Control"
        assert condition_data[1] == "Treatment"
        
        # Test gene data as well
        gene_type_data = mdv.get_column("genes", "gene_type")
        chromosome_data = mdv.get_column("genes", "chromosome")
        
        assert gene_type_data[3] == "ND", "Missing value at index 3 should be 'ND'"
        assert chromosome_data[7] == "ND", "Missing value at index 7 should be 'ND'"
        
    finally:
        if os.path.exists(test_project_dir):
            shutil.rmtree(test_project_dir)


if __name__ == "__main__":
    pytest.main([__file__]) 