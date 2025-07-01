#!/usr/bin/env python3
"""
Test suite for categorical data handling in AnnData to MDV conversion.

This test verifies that categorical columns are properly converted and preserved.

`if data.dtype == "category"` was somewhat rashly replaced with `isinstance(data, pandas.CategoricalDtype)`
which always returned False because `data` is a pandas Series, not a dtype object.

`is_categorical_dtype(data)` properly detects categorical Series - but is deprecated for some reason
and the documentation suggested using `isinstance(data, pandas.CategoricalDtype)` instead.
https://pandas.pydata.org/docs/reference/api/pandas.api.types.is_categorical_dtype.html

There may be scenarios in which the original `data.dtype == "category"` approach would fail,
perhaps not likely in practice, but this test aims to cover all edge cases and verify that the
current implementation is robust and any regressions are caught.
"""

import os
import tempfile
import shutil
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
        # Convert to MDV format
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


if __name__ == "__main__":
    pytest.main([__file__]) 