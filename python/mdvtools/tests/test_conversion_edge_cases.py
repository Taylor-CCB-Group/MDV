#!/usr/bin/env python3
"""
Comprehensive tests for MDV conversion with edge cases and adata.X validation.

This test suite focuses on testing the convert_scanpy_to_mdv function with various
edge cases, ensuring that adata.X is properly handled and validated.
"""

import os
import tempfile
import shutil
import pytest
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from contextlib import contextmanager

from mdvtools.conversions import convert_scanpy_to_mdv
from mdvtools.mdvproject import MDVProject
from mdvtools.spatial.conversion import (
    _compute_table_prefixed_x_umap_and_leiden,
    _concat_spatial_tables,
)
from .mock_anndata import (
    MockAnnDataFactory,
    create_edge_case_anndata,
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


class TestAnnDataXValidation:
    """Test class for validating adata.X structure and handling."""
    
    def test_adata_x_not_none(self):
        """Test that adata.X is not None in created AnnData objects."""
        factory = MockAnnDataFactory(random_seed=42)
        
        # Test all creation methods
        adata_minimal = factory.create_minimal(10, 5)
        adata_realistic = factory.create_realistic(100, 50)
        adata_large = factory.create_large(1000, 200)
        adata_edge = factory.create_edge_cases()
        
        for adata in [adata_minimal, adata_realistic, adata_large, adata_edge]:
            assert adata.X is not None, f"adata.X is None for {type(adata).__name__}"
            assert hasattr(adata.X, 'shape'), f"adata.X has no shape attribute for {type(adata).__name__}"
            assert adata.X.shape == (adata.n_obs, adata.n_vars), f"adata.X shape mismatch for {type(adata).__name__}"
    
    def test_adata_x_type_validation(self):
        """Test that adata.X is either numpy array or scipy sparse matrix."""
        factory = MockAnnDataFactory(random_seed=42)
        
        # Test dense matrix
        adata_dense = factory.create_minimal(10, 5)
        assert isinstance(adata_dense.X, np.ndarray), f"Expected numpy array, got {type(adata_dense.X)}"
        
        # Test sparse matrix
        adata_sparse = factory.create_large(10, 5)  # Large datasets use sparse matrices
        assert isinstance(adata_sparse.X, sp.spmatrix), f"Expected scipy sparse matrix, got {type(adata_sparse.X)}"
    
    def test_adata_x_shape_consistency(self):
        """Test that adata.X shape is consistent with AnnData dimensions."""
        factory = MockAnnDataFactory(random_seed=42)
        
        test_sizes = [(10, 5), (50, 25), (100, 50), (500, 200)]
        
        for n_cells, n_genes in test_sizes:
            adata = factory.create_realistic(n_cells, n_genes)
            expected_shape = (n_cells, n_genes)
            assert adata.X is not None, f"adata.X is None for {type(adata).__name__}"
            actual_shape = adata.X.shape
            
            assert actual_shape == expected_shape, \
                f"Shape mismatch: expected {expected_shape}, got {actual_shape}"
            assert adata.X.shape == (adata.n_obs, adata.n_vars), \
                "adata.X shape doesn't match AnnData dimensions"


class TestConversionWithEdgeCases:
    """Test class for conversion with various edge cases."""
    
    def test_conversion_with_dense_matrix(self):
        """Test conversion with dense numpy array."""
        factory = MockAnnDataFactory(random_seed=42)
        adata = factory.create_minimal(50, 25)
        
        # Ensure it's dense
        assert isinstance(adata.X, np.ndarray)
        
        with temp_mdv_project() as test_dir:
            with suppress_anndata_warnings():
                mdv = convert_scanpy_to_mdv(test_dir, adata, delete_existing=True)
            
            assert isinstance(mdv, MDVProject)
            assert "cells" in mdv.get_datasource_names()
            assert "genes" in mdv.get_datasource_names()
    
    def test_conversion_with_sparse_matrix(self):
        """Test conversion with sparse scipy matrix."""
        factory = MockAnnDataFactory(random_seed=42)
        adata = factory.create_large(100, 50)  # This creates sparse matrix
        
        # Ensure it's sparse
        assert isinstance(adata.X, sp.spmatrix)
        
        with temp_mdv_project() as test_dir:
            with suppress_anndata_warnings():
                mdv = convert_scanpy_to_mdv(test_dir, adata, delete_existing=True)
            
            assert isinstance(mdv, MDVProject)
            assert "cells" in mdv.get_datasource_names()
            assert "genes" in mdv.get_datasource_names()
    
    def test_conversion_with_edge_case_data(self):
        """Test conversion with edge case AnnData."""
        adata = create_edge_case_anndata()
        
        # Validate edge case data
        assert validate_anndata(adata)
        assert adata.X is not None
        
        with temp_mdv_project() as test_dir:
            with suppress_anndata_warnings():
                mdv = convert_scanpy_to_mdv(test_dir, adata, delete_existing=True)
            
            assert isinstance(mdv, MDVProject)
            assert "cells" in mdv.get_datasource_names()
            assert "genes" in mdv.get_datasource_names()
    
    def test_conversion_with_missing_values(self):
        """Test conversion with missing values in metadata."""
        factory = MockAnnDataFactory(random_seed=42)
        adata = factory.create_realistic(100, 50, add_missing=True)
        
        # Check that missing values are present
        summary = get_anndata_summary(adata)
        assert summary['has_missing_obs'] or summary['has_missing_var']
        
        with temp_mdv_project() as test_dir:
            with suppress_anndata_warnings():
                mdv = convert_scanpy_to_mdv(test_dir, adata, delete_existing=True)
            
            assert isinstance(mdv, MDVProject)
    
    def test_conversion_with_large_dataset(self):
        """Test conversion with large dataset to ensure performance."""
        factory = MockAnnDataFactory(random_seed=42)
        adata = factory.create_large(2000, 1000)
        
        # Validate large dataset
        assert validate_anndata(adata)
        assert adata.X is not None
        assert adata.X.shape == (2000, 1000)
        
        with temp_mdv_project() as test_dir:
            with suppress_anndata_warnings():
                mdv = convert_scanpy_to_mdv(test_dir, adata, delete_existing=True)
            
            assert isinstance(mdv, MDVProject)
            assert "cells" in mdv.get_datasource_names()
            assert "genes" in mdv.get_datasource_names()
    
    def test_conversion_with_dimensionality_reductions(self):
        """Test conversion with dimensionality reductions."""
        factory = MockAnnDataFactory(random_seed=42)
        adata = factory.create_realistic(100, 50)
        
        # Check that dimensionality reductions are present
        assert 'X_pca' in adata.obsm
        assert 'X_umap' in adata.obsm
        assert 'PCs' in adata.varm
        
        with temp_mdv_project() as test_dir:
            with suppress_anndata_warnings():
                mdv = convert_scanpy_to_mdv(test_dir, adata, delete_existing=True)
            
            assert isinstance(mdv, MDVProject)
            
            # Check that dimensionality reductions were converted
            cells_metadata = mdv.get_datasource_metadata("cells")
            genes_metadata = mdv.get_datasource_metadata("genes")
            
            # Look for PCA and UMAP columns
            cell_columns = [col['field'] for col in cells_metadata['columns']]
            gene_columns = [col['field'] for col in genes_metadata['columns']]
            
            assert any('pca' in col.lower() for col in cell_columns), "PCA columns not found in cells"
            assert any('umap' in col.lower() for col in cell_columns), "UMAP columns not found in cells"
            assert any('pc' in col.lower() for col in gene_columns), "PC columns not found in genes"
    
    def test_conversion_with_layers(self):
        """Test conversion with expression layers."""
        factory = MockAnnDataFactory(random_seed=42)
        adata = factory.create_realistic(100, 50)
        
        # Check that layers are present
        assert 'counts' in adata.layers
        assert 'log1p' in adata.layers
        assert 'scaled' in adata.layers
        
        with temp_mdv_project() as test_dir:
            with suppress_anndata_warnings():
                mdv = convert_scanpy_to_mdv(test_dir, adata, delete_existing=True)
            
            assert isinstance(mdv, MDVProject)
    
    def test_conversion_with_unstructured_data(self):
        """Test conversion with unstructured data."""
        factory = MockAnnDataFactory(random_seed=42)
        adata = factory.create_realistic(100, 50)
        
        # Check that unstructured data is present
        assert 'neighbors' in adata.uns
        assert 'leiden' in adata.uns
        assert 'rank_genes_groups' in adata.uns
        
        with temp_mdv_project() as test_dir:
            with suppress_anndata_warnings():
                mdv = convert_scanpy_to_mdv(test_dir, adata, delete_existing=True)
            
            assert isinstance(mdv, MDVProject)

    def test_conversion_can_compute_x_umap_and_leiden_as_text(self):
        """Test optional X-based UMAP/Leiden computation during conversion."""
        factory = MockAnnDataFactory(random_seed=42)
        adata = factory.create_minimal(50, 20)
        adata.X = np.asarray(adata.X, dtype=np.float32)
        adata.X[0, 0] = np.nan
        adata.X[1, 1] = np.inf

        assert "X_umap" not in adata.obsm
        assert "leiden" not in adata.obs.columns

        with temp_mdv_project() as test_dir:
            with suppress_anndata_warnings():
                mdv = convert_scanpy_to_mdv(
                    test_dir,
                    adata,
                    delete_existing=True,
                    compute_x_umap=True,
                    leiden_resolution=0.5,
                )

            assert sp.isspmatrix_csc(adata.X)
            assert np.isfinite(adata.X.data).all()

            cells_metadata = mdv.get_datasource_metadata("cells")
            columns = {col["field"]: col for col in cells_metadata["columns"]}

            assert "X_umap_1" in columns
            assert "X_umap_2" in columns
            assert "leiden" in columns
            assert columns["leiden"]["datatype"] in ["text", "text16"]

            leiden_values = mdv.get_column("cells", "leiden")
            assert len(leiden_values) == adata.n_obs
            assert all(isinstance(value, str) for value in leiden_values)

    def test_conversion_can_override_link_name_column_and_datasource_names(self):
        """Test explicit selection of the rows_as_columns link column and datasource names."""
        factory = MockAnnDataFactory(random_seed=42)
        adata = factory.create_minimal(20, 10)
        adata.var["display_name"] = [f"Gene {i}" for i in range(adata.n_vars)]

        with temp_mdv_project() as test_dir:
            with suppress_anndata_warnings():
                mdv = convert_scanpy_to_mdv(
                    test_dir,
                    adata,
                    delete_existing=True,
                    obs_datasource_name="observations",
                    var_datasource_name="features",
                    link_name_column="display_name",
                )

            assert "observations" in mdv.get_datasource_names()
            assert "features" in mdv.get_datasource_names()

            obs_metadata = mdv.get_datasource_metadata("observations")
            link_metadata = obs_metadata["links"]["features"]["rows_as_columns"]
            assert link_metadata["name_column"] == "display_name"

    def test_spatial_single_table_x_umap_is_prefixed_and_text(self):
        """Spatial helper should compute per-table prefixed UMAP/Leiden columns."""
        factory = MockAnnDataFactory(random_seed=42)
        adata = factory.create_minimal(40, 16)

        _compute_table_prefixed_x_umap_and_leiden(
            adata,
            table_prefix="sample_table",
            leiden_resolution=0.6,
        )

        assert "sample_table__X_umap" in adata.obsm
        assert "sample_table__leiden" in adata.obs.columns
        assert "leiden" not in adata.obs.columns
        assert isinstance(adata.obs["sample_table__leiden"].dtype, pd.CategoricalDtype)
        assert all(isinstance(value, str) for value in adata.obs["sample_table__leiden"].astype(str))

        with temp_mdv_project() as test_dir:
            with suppress_anndata_warnings():
                mdv = convert_scanpy_to_mdv(test_dir, adata, delete_existing=True)

            cells_metadata = mdv.get_datasource_metadata("cells")
            columns = {col["field"]: col for col in cells_metadata["columns"]}
            assert "sample_table__X_umap_1" in columns
            assert "sample_table__X_umap_2" in columns
            assert "sample_table__leiden" in columns
            assert columns["sample_table__leiden"]["datatype"] in ["text", "text16"]

    def test_spatial_multi_table_x_umap_stays_per_table_after_merge(self):
        """Merged spatial tables should keep separate per-table UMAP/Leiden columns."""
        factory = MockAnnDataFactory(random_seed=42)
        adata_a = factory.create_minimal(24, 12)
        adata_b = factory.create_minimal(18, 9)

        _compute_table_prefixed_x_umap_and_leiden(
            adata_a,
            table_prefix="table_a",
            leiden_resolution=0.5,
        )
        _compute_table_prefixed_x_umap_and_leiden(
            adata_b,
            table_prefix="table_b",
            leiden_resolution=0.5,
        )
        adata_a.obs["table_name"] = "table_a"
        adata_b.obs["table_name"] = "table_b"

        merged = _concat_spatial_tables([adata_a, adata_b])

        assert "leiden" not in merged.obs.columns
        assert "table_a__leiden" in merged.obs.columns
        assert "table_b__leiden" in merged.obs.columns

        table_a_rows = merged.obs["table_name"] == "table_a"
        table_b_rows = merged.obs["table_name"] == "table_b"

        assert merged.obs.loc[table_a_rows, "table_a__leiden"].notna().all()
        assert merged.obs.loc[table_b_rows, "table_b__leiden"].notna().all()
        assert merged.obs.loc[table_a_rows, "table_b__leiden"].isna().all()
        assert merged.obs.loc[table_b_rows, "table_a__leiden"].isna().all()

        with temp_mdv_project() as test_dir:
            with suppress_anndata_warnings():
                mdv = convert_scanpy_to_mdv(test_dir, merged, delete_existing=True)

            cells_metadata = mdv.get_datasource_metadata("cells")
            columns = {col["field"]: col for col in cells_metadata["columns"]}
            assert "table_a__X_umap_1" in columns
            assert "table_a__X_umap_2" in columns
            assert "table_b__X_umap_1" in columns
            assert "table_b__X_umap_2" in columns
            assert columns["table_a__leiden"]["datatype"] in ["text", "text16"]
            assert columns["table_b__leiden"]["datatype"] in ["text", "text16"]


class TestConversionErrorHandling:
    """Test class for error handling during conversion."""
    
    def test_conversion_with_invalid_adata(self):
        """Test that conversion fails gracefully with invalid AnnData."""
        # Create invalid AnnData with None X
        obs_data = pd.DataFrame({'cell_type': ['A', 'B', 'C']})
        var_data = pd.DataFrame({'gene_type': ['X', 'Y']})
        
        # This should raise an error during conversion
        with pytest.raises(Exception):
            adata = sc.AnnData(X=None, obs=obs_data, var=var_data)
            with temp_mdv_project() as test_dir:
                convert_scanpy_to_mdv(test_dir, adata, delete_existing=True)
    
    def test_conversion_with_shape_mismatch(self):
        """Test that conversion fails with shape mismatch."""
        # Create AnnData with mismatched shapes
        obs_data = pd.DataFrame({'cell_type': ['A', 'B', 'C']})  # 3 cells
        var_data = pd.DataFrame({'gene_type': ['X', 'Y']})      # 2 genes
        X = np.random.random((3, 3))  # 3x3 matrix, should be 3x2
        
        # This should raise an error during AnnData creation
        with pytest.raises(ValueError):
            sc.AnnData(X=X, obs=obs_data, var=var_data)
    
    def test_conversion_with_empty_adata(self):
        """Test conversion with empty AnnData."""
        # Create empty AnnData
        adata = sc.AnnData(X=np.array([]).reshape(0, 0))
        
        # This should fail during conversion
        with pytest.raises(Exception):
            with temp_mdv_project() as test_dir:
                convert_scanpy_to_mdv(test_dir, adata, delete_existing=True)


class TestConversionDataIntegrity:
    """Test class for ensuring data integrity during conversion."""
    
    def test_conversion_preserves_data_structure(self):
        """Test that conversion preserves the basic data structure."""
        factory = MockAnnDataFactory(random_seed=42)
        adata = factory.create_realistic(100, 50)
        
        # Get original data summary
        original_summary = get_anndata_summary(adata)
        
        with temp_mdv_project() as test_dir:
            with suppress_anndata_warnings():
                mdv = convert_scanpy_to_mdv(test_dir, adata, delete_existing=True)
            
            # Verify basic structure is preserved
            assert mdv.get_datasource_metadata("cells") is not None
            assert mdv.get_datasource_metadata("genes") is not None
            
            # Check that number of columns matches expectations
            cells_metadata = mdv.get_datasource_metadata("cells")
            genes_metadata = mdv.get_datasource_metadata("genes")
            
            # Should have at least the original obs/var columns plus some computed ones
            assert len(cells_metadata['columns']) >= len(original_summary['obs_columns'])
            assert len(genes_metadata['columns']) >= len(original_summary['var_columns'])
    
    def test_conversion_preserves_categorical_data(self):
        """Test that categorical data is preserved during conversion."""
        factory = MockAnnDataFactory(random_seed=42)
        adata = factory.create_realistic(100, 50)
        
        # Check original categorical columns
        original_categorical_obs = [col for col in adata.obs.columns 
                                  if hasattr(adata.obs[col], 'cat')]
        original_categorical_var = [col for col in adata.var.columns 
                                  if hasattr(adata.var[col], 'cat')]
        
        with temp_mdv_project() as test_dir:
            with suppress_anndata_warnings():
                mdv = convert_scanpy_to_mdv(test_dir, adata, delete_existing=True)
            
            # Check that categorical columns are converted to text type
            cells_metadata = mdv.get_datasource_metadata("cells")
            genes_metadata = mdv.get_datasource_metadata("genes")
            
            for col_name in original_categorical_obs:
                col_meta = next((col for col in cells_metadata['columns'] 
                               if col['field'] == col_name), None)
                assert col_meta is not None, f"Categorical column {col_name} not found in cells"
                assert col_meta['datatype'] in ['text', 'text16'], \
                    f"Categorical column {col_name} should be text type, got {col_meta['datatype']}"
            
            for col_name in original_categorical_var:
                col_meta = next((col for col in genes_metadata['columns'] 
                               if col['field'] == col_name), None)
                assert col_meta is not None, f"Categorical column {col_name} not found in genes"
                assert col_meta['datatype'] in ['text', 'text16'], \
                    f"Categorical column {col_name} should be text type, got {col_meta['datatype']}"


if __name__ == "__main__":
    pytest.main([__file__]) 
