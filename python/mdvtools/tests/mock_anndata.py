#!/usr/bin/env python3
"""
Mock AnnData generation module for testing and stress testing.

This module provides utilities to create realistic AnnData objects with various
configurations, data types, and edge cases for comprehensive testing of the
MDV conversion pipeline.
"""

import numpy as np
import pandas as pd
import scanpy as sc
import scipy
import scipy.sparse
from typing import Dict, List, Optional, Union, Any
import warnings
from contextlib import contextmanager


@contextmanager
def suppress_anndata_warnings():
    """Context manager to suppress expected AnnData warnings."""
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Transforming to str index", category=UserWarning)
        yield


def chunked_log1p_normalization(sparse_matrix, chunk_size=1000):
    """Perform log1p normalization in chunks to avoid dense matrices."""
    if not scipy.sparse.issparse(sparse_matrix):
        # For dense matrices, just apply log1p directly
        return np.log1p(sparse_matrix)
    n_cells, n_genes = sparse_matrix.shape
    if n_cells <= chunk_size:
        # preserve sparsity
        result = sparse_matrix.copy()
        result.data = np.log1p(result.data)
        return result
    # For large matrices, process in chunks and vstack
    chunks = []
    for chunk_start in range(0, n_cells, chunk_size):
        chunk_end = min(chunk_start + chunk_size, n_cells)
        chunk = sparse_matrix[chunk_start:chunk_end, :].toarray()
        chunk_log = np.log1p(chunk)
        chunk_sparse = scipy.sparse.csr_matrix(chunk_log)
        chunks.append(chunk_sparse)
    return scipy.sparse.vstack(chunks, format='csr')


def chunked_zscore_normalization(sparse_matrix, chunk_size=1000):
    """Perform z-score normalization in chunks to avoid dense matrices."""
    if not scipy.sparse.issparse(sparse_matrix):
        X_dense = np.asarray(sparse_matrix)
        return (X_dense - X_dense.mean(axis=0)) / (X_dense.std(axis=0) + 1e-8)
    n_cells, n_genes = sparse_matrix.shape
    if n_cells <= chunk_size:
        X_dense = sparse_matrix.toarray()
        return (X_dense - X_dense.mean(axis=0)) / (X_dense.std(axis=0) + 1e-8)
    # Compute mean and std across all rows
    total_sum = np.zeros(n_genes)
    total_sum_sq = np.zeros(n_genes)
    for chunk_start in range(0, n_cells, chunk_size):
        chunk_end = min(chunk_start + chunk_size, n_cells)
        chunk = sparse_matrix[chunk_start:chunk_end, :].toarray()
        total_sum += chunk.sum(axis=0)
        total_sum_sq += (chunk ** 2).sum(axis=0)
    mean_vals = total_sum / n_cells
    var_vals = (total_sum_sq / n_cells) - (mean_vals ** 2)
    std_vals = np.sqrt(var_vals + 1e-8)
    # Now apply normalization in chunks and vstack
    chunks = []
    for chunk_start in range(0, n_cells, chunk_size):
        chunk_end = min(chunk_start + chunk_size, n_cells)
        chunk = sparse_matrix[chunk_start:chunk_end, :].toarray()
        chunk_normalized = (chunk - mean_vals) / std_vals
        chunk_sparse = scipy.sparse.csr_matrix(chunk_normalized)
        chunks.append(chunk_sparse)
    return scipy.sparse.vstack(chunks, format='csr')


def estimate_memory_usage(n_cells, n_genes, sparse=True):
    """Estimate memory usage for a dataset.
    
    Args:
        n_cells: Number of cells
        n_genes: Number of genes
        sparse: Whether the matrix is sparse
        
    Returns:
        Estimated memory usage in MB
    """
    if sparse:
        # Estimate 10% sparsity for typical single-cell data
        sparsity = 0.1
        nnz = int(n_cells * n_genes * sparsity)
        # 8 bytes per value + 4 bytes per index + 4 bytes per indptr
        memory_bytes = nnz * 16 + n_cells * 4
    else:
        # Dense matrix: 8 bytes per element
        memory_bytes = n_cells * n_genes * 8
    
    return memory_bytes / (1024 * 1024)  # Convert to MB


class MockAnnDataFactory:
    """Factory class for creating mock AnnData objects with various configurations."""
    
    def __init__(self, random_seed: Optional[int] = None):
        """Initialize the factory with optional random seed."""
        if random_seed is not None:
            np.random.seed(random_seed)
    
    def create_minimal(self, n_cells: int = 10, n_genes: int = 5, 
                      add_missing: bool = False) -> sc.AnnData:
        """Create a minimal AnnData object for basic testing."""
        return self._create_anndata(
            n_cells=n_cells, 
            n_genes=n_genes, 
            add_missing=add_missing,
            add_dim_reductions=False,
            add_layers=False,
            add_uns=False
        )
    
    def create_realistic(self, n_cells: int = 1000, n_genes: int = 2000,
                        add_missing: bool = True) -> sc.AnnData:
        """Create a realistic AnnData object with typical single-cell data features."""
        return self._create_anndata(
            n_cells=n_cells,
            n_genes=n_genes,
            add_missing=add_missing,
            add_dim_reductions=True,
            add_layers=True,
            add_uns=True
        )
    
    def create_large(self, n_cells: int = 10000, n_genes: int = 5000,
                    add_missing: bool = True) -> sc.AnnData:
        """Create a large AnnData object for stress testing."""
        return self._create_anndata(
            n_cells=n_cells,
            n_genes=n_genes,
            add_missing=add_missing,
            add_dim_reductions=True,
            add_layers=True,
            add_uns=True,
            sparse_matrix=True
        )
    
    def create_memory_efficient_large(self, n_cells: int = 10000, n_genes: int = 5000,
                                    add_missing: bool = True) -> sc.AnnData:
        """Create a large AnnData object optimized for memory efficiency.
        
        This method creates large datasets without dense layers to avoid
        excessive memory consumption during stress testing.
        """
        return self._create_anndata(
            n_cells=n_cells,
            n_genes=n_genes,
            add_missing=add_missing,
            add_dim_reductions=False,  # Skip dense dimensionality reductions
            add_layers=False,          # Skip dense layers
            add_uns=True,
            sparse_matrix=True
        )
    
    def create_massive_dataset(self, n_cells: int = 100000, n_genes: int = 10000,
                             add_missing: bool = True) -> sc.AnnData:
        """Create a massive dataset (100k+ cells) for extreme stress testing.
        
        This method uses chunked operations and memory-efficient approaches
        to handle datasets that would otherwise cause memory issues.
        """
        print(f"Creating massive dataset: {n_cells:,} cells x {n_genes:,} genes")
        print(f"Estimated memory usage: {estimate_memory_usage(n_cells, n_genes, sparse=True):.1f}MB (sparse)")
        
        return self._create_anndata(
            n_cells=n_cells,
            n_genes=n_genes,
            add_missing=add_missing,
            add_dim_reductions=False,  # Skip dense dimensionality reductions
            add_layers=True,           # Use chunked layers
            add_uns=True,
            sparse_matrix=True,
            use_chunked_layers=True    # Enable chunked layer processing
        )
    
    def create_edge_cases(self) -> sc.AnnData:
        """Create an AnnData object with various edge cases and problematic data."""
        return self._create_edge_case_anndata()
    
    def create_with_specific_features(self, 
                                    cell_types: Optional[List[str]] = None,
                                    conditions: Optional[List[str]] = None,
                                    gene_types: Optional[List[str]] = None,
                                    n_cells: int = 100,
                                    n_genes: int = 200,
                                    **kwargs) -> sc.AnnData:
        """Create AnnData with specific categorical features."""
        return self._create_anndata(
            n_cells=n_cells,
            n_genes=n_genes,
            cell_types=cell_types,
            conditions=conditions,
            gene_types=gene_types,
            **kwargs
        )
    
    def _create_anndata(self, n_cells: int, n_genes: int, 
                       add_missing: bool = False,
                       add_dim_reductions: bool = False,
                       add_layers: bool = False,
                       add_uns: bool = False,
                       sparse_matrix: bool = False,
                       cell_types: Optional[List[str]] = None,
                       conditions: Optional[List[str]] = None,
                       gene_types: Optional[List[str]] = None,
                       use_chunked_layers: bool = False) -> sc.AnnData:
        """Internal method to create AnnData with specified features."""
        
        # Default categorical values
        if cell_types is None:
            cell_types = ['T-cell', 'B-cell', 'NK-cell', 'Monocyte', 'Dendritic']
        if conditions is None:
            conditions = ['Control', 'Treatment', 'Disease']
        if gene_types is None:
            gene_types = ['protein_coding', 'lncRNA', 'miRNA', 'pseudogene', 'rRNA']
        
        # Create cell metadata (obs)
        obs_data = self._create_obs_data(
            n_cells, cell_types, conditions, add_missing
        )
        
        # Create gene metadata (var)
        var_data = self._create_var_data(
            n_genes, gene_types, add_missing
        )
        
        # Create expression matrix
        X = self._create_expression_matrix(n_cells, n_genes, sparse_matrix)
        
        # Create AnnData object
        adata = sc.AnnData(X=X, obs=obs_data, var=var_data)
        
        # Add dimensionality reductions
        if add_dim_reductions:
            self._add_dimension_reductions(adata)
        
        # Add layers
        if add_layers:
            self._add_layers(adata, use_chunked_layers)
        
        # Add unstructured data
        if add_uns:
            self._add_unstructured_data(adata)
        
        return adata
    
    def _create_obs_data(self, n_cells: int, cell_types: List[str], 
                        conditions: List[str], add_missing: bool) -> pd.DataFrame:
        """Create cell metadata DataFrame."""
        # Generate probability arrays that match the number of categories
        cell_type_probs = [1.0 / len(cell_types)] * len(cell_types)
        condition_probs = [1.0 / len(conditions)] * len(conditions)
        
        obs_data = pd.DataFrame({
            'cell_type': pd.Categorical(
                np.random.choice(cell_types, n_cells, p=cell_type_probs)
            ),
            'condition': pd.Categorical(
                np.random.choice(conditions, n_cells, p=condition_probs)
            ),
            'quality_score': np.random.normal(0, 1, n_cells),
            'total_counts': np.random.exponential(1000, n_cells),
            'n_genes_by_counts': np.random.poisson(2000, n_cells),
            'pct_counts_mt': np.random.beta(2, 20, n_cells) * 10,
            'is_high_quality': pd.Series(
                np.random.choice([True, False], n_cells, p=[0.8, 0.2]), 
                dtype='object'
            ),
            'is_doublet': pd.Series(
                np.random.choice([True, False], n_cells, p=[0.1, 0.9]), 
                dtype='object'
            ),
            'patient_id': pd.Categorical(
                [f'P{i:03d}' for i in np.random.randint(1, 21, n_cells)]
            ),
            'batch': pd.Categorical(
                [f'batch_{i}' for i in np.random.randint(1, 6, n_cells)]
            )
        })
        
        # Add missing values if requested
        if add_missing:
            missing_indices = np.random.choice(n_cells, size=n_cells//10, replace=False)
            for idx in missing_indices:
                col = np.random.choice(['cell_type', 'condition', 'is_high_quality', 'is_doublet'])
                obs_data.loc[idx, col] = np.nan
        
        return obs_data
    
    def _create_var_data(self, n_genes: int, gene_types: List[str], 
                        add_missing: bool) -> pd.DataFrame:
        """Create gene metadata DataFrame."""
        # Generate probability arrays that match the number of categories
        gene_type_probs = [1.0 / len(gene_types)] * len(gene_types)
        
        var_data = pd.DataFrame({
            'gene_type': pd.Categorical(
                np.random.choice(gene_types, n_genes, p=gene_type_probs)
            ),
            'chromosome': pd.Categorical(
                [f'chr{i}' for i in np.random.randint(1, 23, n_genes)]
            ),
            'mean_expression': np.random.exponential(1, n_genes),
            'highly_variable': np.random.choice([True, False], n_genes, p=[0.2, 0.8]),
            'mt': pd.Series(
                [name.startswith('MT-') for name in [f'GENE_{i:05d}' for i in range(n_genes)]],
                dtype='object'
            ),
            'ribosomal': pd.Series(
                [name.startswith('RPS') or name.startswith('RPL') 
                 for name in [f'GENE_{i:05d}' for i in range(n_genes)]],
                dtype='object'
            ),
            'name': [f'GENE_{i:05d}' for i in range(n_genes)]
        })
        
        # Add missing values if requested
        if add_missing:
            missing_indices = np.random.choice(n_genes, size=n_genes//10, replace=False)
            for idx in missing_indices:
                col = np.random.choice(['gene_type', 'highly_variable', 'mt', 'ribosomal'])
                var_data.loc[idx, col] = np.nan
        
        return var_data
    
    def _create_expression_matrix(self, n_cells: int, n_genes: int, 
                                sparse: bool = False) -> Union[np.ndarray, scipy.sparse.spmatrix]:
        """Create expression matrix with realistic single-cell data patterns."""
        if sparse:
            # Create sparse matrix with realistic sparsity
            sparsity = 0.9  # 90% zeros
            nnz = int(n_cells * n_genes * (1 - sparsity))
            
            # Generate random indices
            cell_indices = np.random.randint(0, n_cells, nnz)
            gene_indices = np.random.randint(0, n_genes, nnz)
            
            # Generate expression values (negative binomial distribution)
            values = np.random.negative_binomial(5, 0.3, nnz)
            
            # Create sparse matrix
            import scipy.sparse as sp
            X = sp.csr_matrix((values, (cell_indices, gene_indices)), 
                            shape=(n_cells, n_genes))
        else:
            # Create dense matrix
            X = np.random.negative_binomial(5, 0.3, (n_cells, n_genes))
        
        return X
    
    def _add_dimension_reductions(self, adata: sc.AnnData):
        """Add dimensionality reductions to the AnnData object."""
        n_cells = adata.n_obs
        n_genes = adata.n_vars
        
        # PCA for cells
        n_pcs = min(50, n_cells, n_genes)
        pca_cells = np.random.normal(0, 1, (n_cells, n_pcs))
        adata.obsm['X_pca'] = pca_cells
        
        # UMAP for cells
        umap_cells = np.random.normal(0, 1, (n_cells, 2))
        adata.obsm['X_umap'] = umap_cells
        
        # t-SNE for cells
        tsne_cells = np.random.normal(0, 1, (n_cells, 2))
        adata.obsm['X_tsne'] = tsne_cells
        
        # PCA for genes
        pca_genes = np.random.normal(0, 1, (n_genes, n_pcs))
        adata.varm['PCs'] = pca_genes
    
    def _add_layers(self, adata: sc.AnnData, use_chunked_layers: bool):
        """Add expression layers to the AnnData object."""
        # Ensure adata.X is valid
        assert adata.X is not None, "adata.X cannot be None"
        assert hasattr(adata.X, 'shape'), "adata.X must have a shape attribute"
        assert adata.X.shape == (adata.n_obs, adata.n_vars), f"adata.X shape {adata.X.shape} doesn't match AnnData dimensions ({adata.n_obs}, {adata.n_vars})"
        
        # Raw counts - preserve sparsity if input is sparse
        if scipy.sparse.issparse(adata.X):
            # For sparse matrices, preserve sparsity by copying the sparse structure
            adata.layers['counts'] = adata.X.copy()  # type: ignore
        else:
            # For dense arrays, use numpy copy
            adata.layers['counts'] = np.array(adata.X, copy=True)
        
        # Log-normalized data - use chunked normalization if requested
        if use_chunked_layers:
            log_data = chunked_log1p_normalization(adata.X)
        else:
            if scipy.sparse.issparse(adata.X):
                # preserve sparsity
                log_data = scipy.sparse.csc_matrix(adata.X.copy()) # type: ignore
                log_data.data = np.log1p(log_data.data)
            else:
                log_data = np.log1p(np.asarray(adata.X))
        
        adata.layers['log1p'] = log_data # type: ignore
        
        # Scaled data - use chunked normalization for large datasets
        if use_chunked_layers and scipy.sparse.issparse(log_data):
            # Use chunked normalization directly on sparse matrix
            scaled_data = chunked_zscore_normalization(log_data)
        else:
            # For smaller datasets or when chunking is disabled, use traditional approach
            if scipy.sparse.issparse(log_data):
                log_dense = log_data.toarray()  # type: ignore
            else:
                log_dense = np.asarray(log_data)
            
            if use_chunked_layers:
                scaled_data = chunked_zscore_normalization(log_dense)
            else:
                scaled_data = (log_dense - log_dense.mean(axis=0)) / (log_dense.std(axis=0) + 1e-8) # type: ignore
        
        adata.layers['scaled'] = scaled_data # type: ignore
    
    def _add_unstructured_data(self, adata: sc.AnnData):
        """Add unstructured data to the AnnData object."""
        adata.uns['neighbors'] = {
            'params': {'n_neighbors': 15, 'metric': 'euclidean'},
            'connectivities_key': 'connectivities',
            'distances_key': 'distances'
        }
        
        adata.uns['leiden'] = {
            'params': {'resolution': 0.5, 'random_state': 0},
            'connectivities_key': 'connectivities'
        }
        
        adata.uns['rank_genes_groups'] = {
            'params': {'groupby': 'leiden', 'method': 't-test'},
            'names': [['gene1', 'gene2', 'gene3'] for _ in range(5)],
            'scores': [[1.5, 1.2, 1.0] for _ in range(5)],
            'logfoldchanges': [[0.8, 0.6, 0.4] for _ in range(5)]
        }
    
    def _create_edge_case_anndata(self) -> sc.AnnData:
        """Create AnnData with various edge cases and problematic data."""
        n_cells, n_genes = 50, 100
        
        # Create obs with edge cases - ensure all arrays have the same length
        special_chars_pattern = ['test\nnewline', 'test\ttab', 'test"quote', 'test\'apos']
        unicode_pattern = ['αβγδε', '🎉🎊🎈', '中文测试']
        mixed_types_pattern = [1, 'string', True, 3.14, None]
        inf_values_pattern = [float('inf'), float('-inf')]
        boolean_pattern = [True, False, np.nan]
        
        # Repeat patterns to match n_cells
        special_chars = (special_chars_pattern * ((n_cells // len(special_chars_pattern)) + 1))[:n_cells]
        unicode_chars = (unicode_pattern * ((n_cells // len(unicode_pattern)) + 1))[:n_cells]
        mixed_types = (mixed_types_pattern * ((n_cells // len(mixed_types_pattern)) + 1))[:n_cells]
        inf_values = (inf_values_pattern * ((n_cells // len(inf_values_pattern)) + 1))[:n_cells]
        boolean_values = (boolean_pattern * ((n_cells // len(boolean_pattern)) + 1))[:n_cells]
        
        obs_data = pd.DataFrame({
            'empty_string': [''] * n_cells,
            'very_long_string': ['A' * 1000] * n_cells,
            'special_chars': special_chars,
            'unicode': unicode_chars,
            'mixed_types': mixed_types,
            'all_nan': [np.nan] * n_cells,
            'inf_values': inf_values,
            'zero_variance': [42] * n_cells,
            'boolean_with_nan': pd.Series(boolean_values, dtype='object')
        })
        
        # Create var with edge cases - ensure all arrays have the same length
        mixed_dtypes_pattern = [1, 'string', True, 3.14]
        mixed_dtypes = (mixed_dtypes_pattern * ((n_genes // len(mixed_dtypes_pattern)) + 1))[:n_genes]
        
        var_data = pd.DataFrame({
            'empty_categories': pd.Categorical([''] * n_genes),
            'single_category': pd.Categorical(['same'] * n_genes),
            'many_categories': pd.Categorical([f'cat_{i}' for i in range(n_genes)]),
            'mixed_dtypes': mixed_dtypes,
            'all_inf': [float('inf')] * n_genes,
            'all_nan': [np.nan] * n_genes
        })
        
        # Create problematic expression matrix
        X = np.random.normal(0, 1, (n_cells, n_genes))
        # Add some problematic values
        X[0, 0] = float('inf')
        X[0, 1] = float('-inf')
        X[0, 2] = np.nan
        
        return sc.AnnData(X=X, obs=obs_data, var=var_data)


# Convenience functions for backward compatibility
def create_minimal_anndata(n_cells: int = 10, n_genes: int = 5, 
                          add_missing: bool = False) -> sc.AnnData:
    """Create minimal AnnData object for testing (backward compatibility).
    
    This function maintains exact backward compatibility with the original implementation
    to ensure existing tests continue to pass.
    """
    # Cell metadata - use the original simple pattern
    obs_data = pd.DataFrame({
        'cell_type': pd.Categorical(['T-cell', 'B-cell'] * (n_cells // 2)),
        'condition': pd.Categorical(['Control', 'Treatment'] * (n_cells // 2)),
        'quality_score': np.random.normal(0, 1, n_cells)
    })
    
    # Gene metadata - use the original simple pattern  
    var_data = pd.DataFrame({
        'gene_type': pd.Categorical(['protein_coding', 'lncRNA'] * (n_genes // 2) + ['miRNA']),
        'chromosome': pd.Categorical(['chr1', 'chr2'] * (n_genes // 2) + ['chrX']),
        'mean_expression': np.random.exponential(1, n_genes)
    })
    
    # Add missing values if requested - use the original deterministic placement
    if add_missing:
        obs_data.loc[2, 'cell_type'] = np.nan
        obs_data.loc[5, 'condition'] = np.nan
        var_data.loc[3, 'gene_type'] = np.nan
    
    # Expression matrix - use the original distribution
    X = np.random.negative_binomial(5, 0.3, (n_cells, n_genes))
    
    return sc.AnnData(X=X, obs=obs_data, var=var_data)


def create_realistic_anndata(n_cells: int = 1000, n_genes: int = 2000,
                           add_missing: bool = True) -> sc.AnnData:
    """Create realistic AnnData object with typical single-cell data features."""
    factory = MockAnnDataFactory()
    return factory.create_realistic(n_cells, n_genes, add_missing)


def create_large_anndata(n_cells: int = 10000, n_genes: int = 5000,
                        add_missing: bool = True) -> sc.AnnData:
    """Create large AnnData object for stress testing."""
    factory = MockAnnDataFactory()
    return factory.create_large(n_cells, n_genes, add_missing)


def create_memory_efficient_large_anndata(n_cells: int = 10000, n_genes: int = 5000,
                                         add_missing: bool = True) -> sc.AnnData:
    """Create memory-efficient large AnnData object for stress testing.
    
    This function creates large datasets without dense layers to avoid
    excessive memory consumption during stress testing.
    """
    factory = MockAnnDataFactory()
    return factory.create_memory_efficient_large(n_cells, n_genes, add_missing)


def create_massive_anndata(n_cells: int = 100000, n_genes: int = 10000,
                          add_missing: bool = True) -> sc.AnnData:
    """Create massive AnnData object for extreme stress testing.
    
    This function creates datasets with 100k+ cells using chunked operations
    to handle memory efficiently. Suitable for testing with real-world scale data.
    """
    factory = MockAnnDataFactory()
    return factory.create_massive_dataset(n_cells, n_genes, add_missing)


def create_edge_case_anndata() -> sc.AnnData:
    """Create AnnData object with various edge cases and problematic data."""
    factory = MockAnnDataFactory()
    return factory.create_edge_cases()


# Utility functions for testing
def get_anndata_summary(adata: sc.AnnData) -> Dict[str, Any]:
    """Get a summary of AnnData object properties for testing."""
    # Ensure adata.X is valid
    assert adata.X is not None, "adata.X cannot be None"
    
    return {
        'n_cells': adata.n_obs,
        'n_genes': adata.n_vars,
        'obs_columns': list(adata.obs.columns),
        'var_columns': list(adata.var.columns),
        'obsm_keys': list(adata.obsm.keys()),
        'varm_keys': list(adata.varm.keys()),
        'layers_keys': list(adata.layers.keys()),
        'uns_keys': list(adata.uns.keys()),
        'sparse': hasattr(adata.X, 'toarray'),
        'has_missing_obs': bool(adata.obs.isnull().values.any()),
        'has_missing_var': bool(adata.var.isnull().values.any()),
        'categorical_obs': [col for col in adata.obs.columns 
                          if hasattr(adata.obs[col], 'cat')],
        'categorical_var': [col for col in adata.var.columns 
                          if hasattr(adata.var[col], 'cat')]
    }


def validate_anndata(adata: sc.AnnData) -> bool:
    """Validate that AnnData object has expected structure."""
    try:
        # Basic structure
        assert hasattr(adata, 'obs') and hasattr(adata, 'var') and hasattr(adata, 'X')
        assert adata.n_obs > 0 and adata.n_vars > 0
        
        # Check that obs and var have the right number of rows
        assert len(adata.obs) == adata.n_obs
        assert len(adata.var) == adata.n_vars
        
        # Check that X is not None and has the right shape
        assert adata.X is not None, "adata.X cannot be None"
        assert hasattr(adata.X, 'shape'), "adata.X must have a shape attribute"
        assert adata.X.shape == (adata.n_obs, adata.n_vars), f"adata.X shape {adata.X.shape} doesn't match AnnData dimensions ({adata.n_obs}, {adata.n_vars})"
        
        return True
    except AssertionError:
        return False 