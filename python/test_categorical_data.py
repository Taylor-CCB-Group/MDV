#!/usr/bin/env python3
"""
Test script for categorical data handling in AnnData to MDV conversion.
This test verifies that categorical columns are properly converted and preserved.
"""

import sys
import os
import tempfile
import shutil
import pandas as pd
import numpy as np
import scanpy as sc

# Add the python directory to the path so we can import mdvtools
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

from mdvtools.conversions import convert_scanpy_to_mdv
from mdvtools.mdvproject import MDVProject


def create_test_anndata_with_categorical():
    """Create a synthetic AnnData object with categorical data for testing."""
    # Create synthetic data
    n_cells = 100
    n_genes = 50
    
    # Create cell metadata with categorical columns
    cell_types = ['T-cell', 'B-cell', 'NK-cell', 'Monocyte']
    conditions = ['Control', 'Treatment']
    batches = ['Batch1', 'Batch2', 'Batch3']
    
    obs_data = pd.DataFrame({
        'cell_type': pd.Categorical(np.random.choice(cell_types, n_cells)),
        'condition': pd.Categorical(np.random.choice(conditions, n_cells)),
        'batch': pd.Categorical(np.random.choice(batches, n_cells)),
        'quality_score': np.random.normal(0, 1, n_cells),  # numeric column
        'is_high_quality': pd.Categorical(np.random.choice(['Yes', 'No'], n_cells))
    })
    
    # Create gene metadata with categorical columns
    gene_types = ['protein_coding', 'lncRNA', 'pseudogene', 'miRNA']
    chromosomes = ['chr1', 'chr2', 'chrX', 'chrY']
    
    var_data = pd.DataFrame({
        'gene_type': pd.Categorical(np.random.choice(gene_types, n_genes)),
        'chromosome': pd.Categorical(np.random.choice(chromosomes, n_genes)),
        'is_mitochondrial': pd.Categorical(np.random.choice(['Yes', 'No'], n_genes)),
        'mean_expression': np.random.exponential(1, n_genes)  # numeric column
    })
    
    # Create expression matrix
    X = np.random.negative_binomial(5, 0.3, (n_cells, n_genes))
    
    # Create AnnData object
    adata = sc.AnnData(X=X, obs=obs_data, var=var_data)
    
    return adata


def test_categorical_data_conversion():
    """Test that categorical data is properly converted from AnnData to MDV format."""
    print("Creating test AnnData with categorical data...")
    
    # Create test AnnData with categorical data
    adata = create_test_anndata_with_categorical()
    
    # Create temporary directory for testing
    test_project_dir = tempfile.mkdtemp()
    
    try:
        print(f"Converting to MDV format in {test_project_dir}...")
        
        # Convert to MDV format
        mdv = convert_scanpy_to_mdv(test_project_dir, adata, delete_existing=True)
        
        # Basic assertions
        assert isinstance(mdv, MDVProject), "MDV project should be created"
        print("✓ MDV project created successfully")
        
        # Check if datasources were created
        datasource_names = mdv.get_datasource_names()
        assert "cells" in datasource_names, "Cells datasource not found"
        assert "genes" in datasource_names, "Genes datasource not found"
        print("✓ Cells and genes datasources created")
        
        # Get datasource metadata
        cell_datasource = mdv.get_datasource_metadata("cells")
        gene_datasource = mdv.get_datasource_metadata("genes")
        
        # Verify categorical columns in cells datasource
        cell_columns = {col['field']: col for col in cell_datasource['columns']}
        
        # Check that categorical columns are properly handled
        categorical_cell_cols = ['cell_type', 'condition', 'batch', 'is_high_quality']
        for col_name in categorical_cell_cols:
            assert col_name in cell_columns, f"Categorical column {col_name} not found in cells datasource"
            col_meta = cell_columns[col_name]
            assert col_meta['datatype'] in ['text', 'text16'], f"Column {col_name} should be text type, got {col_meta['datatype']}"
            assert 'values' in col_meta, f"Column {col_name} should have values list"
            assert len(col_meta['values']) > 0, f"Column {col_name} should have non-empty values list"
            print(f"✓ Categorical column '{col_name}' properly handled in cells datasource")
        
        # Check numeric columns are handled correctly
        numeric_cell_cols = ['quality_score']
        for col_name in numeric_cell_cols:
            assert col_name in cell_columns, f"Numeric column {col_name} not found in cells datasource"
            col_meta = cell_columns[col_name]
            assert col_meta['datatype'] in ['double', 'integer'], f"Column {col_name} should be numeric type, got {col_meta['datatype']}"
            print(f"✓ Numeric column '{col_name}' properly handled in cells datasource")
        
        # Verify categorical columns in genes datasource
        gene_columns = {col['field']: col for col in gene_datasource['columns']}
        
        categorical_gene_cols = ['gene_type', 'chromosome', 'is_mitochondrial']
        for col_name in categorical_gene_cols:
            assert col_name in gene_columns, f"Categorical column {col_name} not found in genes datasource"
            col_meta = gene_columns[col_name]
            assert col_meta['datatype'] in ['text', 'text16'], f"Column {col_name} should be text type, got {col_meta['datatype']}"
            assert 'values' in col_meta, f"Column {col_name} should have values list"
            assert len(col_meta['values']) > 0, f"Column {col_name} should have non-empty values list"
            print(f"✓ Categorical column '{col_name}' properly handled in genes datasource")
        
        # Check numeric columns in genes datasource
        numeric_gene_cols = ['mean_expression']
        for col_name in numeric_gene_cols:
            assert col_name in gene_columns, f"Numeric column {col_name} not found in genes datasource"
            col_meta = gene_columns[col_name]
            assert col_meta['datatype'] in ['double', 'integer'], f"Column {col_name} should be numeric type, got {col_meta['datatype']}"
            print(f"✓ Numeric column '{col_name}' properly handled in genes datasource")
        
        # Test that we can retrieve the data correctly
        for col_name in categorical_cell_cols:
            retrieved_data = mdv.get_column("cells", col_name)
            assert len(retrieved_data) == 100, f"Retrieved data for {col_name} should have 100 rows"
            # Check that all values are strings (categorical data should be converted to text)
            assert all(isinstance(val, str) for val in retrieved_data), f"All values in {col_name} should be strings"
            print(f"✓ Data retrieval works for categorical column '{col_name}'")
        
        # Test that unique values are preserved
        for col_name in categorical_cell_cols:
            retrieved_data = mdv.get_column("cells", col_name)
            original_unique = set(adata.obs[col_name].cat.categories)
            retrieved_unique = set(retrieved_data)
            # Note: ND might be added for missing values, so we check that original values are preserved
            assert original_unique.issubset(retrieved_unique), f"Original categorical values for {col_name} should be preserved"
            print(f"✓ Original categorical values preserved for '{col_name}'")
        
        print("\n🎉 All categorical data tests passed!")
        return True
        
    except Exception as e:
        print(f"\n❌ Test failed with error: {e}")
        import traceback
        traceback.print_exc()
        return False
        
    finally:
        # Clean up
        if os.path.exists(test_project_dir):
            shutil.rmtree(test_project_dir)
            print(f"Cleaned up test directory: {test_project_dir}")


if __name__ == "__main__":
    success = test_categorical_data_conversion()
    sys.exit(0 if success else 1) 