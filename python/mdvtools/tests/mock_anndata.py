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
                    add_missing: bool = True, density: float = 0.1) -> sc.AnnData:
        """Create a large AnnData object for stress testing."""
        return self._create_anndata(
            n_cells=n_cells,
            n_genes=n_genes,
            add_missing=add_missing,
            add_dim_reductions=True,
            add_layers=True,
            add_uns=True,
            sparse_matrix=True,
            density=density
        )
    
    def create_memory_efficient_large(self, n_cells: int = 10000, n_genes: int = 5000,
                                    add_missing: bool = True, density: float = 0.1) -> sc.AnnData:
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
            add_uns=False,             # Skip unstructured data for memory efficiency
            sparse_matrix=True,
            density=density,
            minimal_metadata=True      # Use minimal metadata for large datasets
        )
    
    def create_massive_dataset(self, n_cells: int = 100000, n_genes: int = 10000,
                             add_missing: bool = True, density: float = 0.1,
                             chunk_size: int = 10000, mode: str = 'realistic') -> sc.AnnData:
        """Create a massive dataset (100k+ cells) for extreme stress testing.
        
        This method uses chunked operations and memory-efficient approaches
        to handle datasets that would otherwise cause memory issues.
        
        Args:
            n_cells: Number of cells
            n_genes: Number of genes
            add_missing: Whether to add missing values
            density: Density of non-zero elements
            chunk_size: Size of chunks for matrix generation
            mode: Generation mode - 'realistic', 'fast', or 'skeleton'
        """
        print(f"Creating massive dataset: {n_cells:,} cells x {n_genes:,} genes")
        print(f"Estimated memory usage: {estimate_memory_usage(n_cells, n_genes, sparse=True):.1f}MB (sparse)")
        print(f"Mode: {mode}, Chunk size: {chunk_size:,}")
        
        return self._create_anndata(
            n_cells=n_cells,
            n_genes=n_genes,
            add_missing=add_missing,
            add_dim_reductions=False,  # Skip dense dimensionality reductions
            add_layers=False,          # Skip layers for maximum memory efficiency
            add_uns=False,             # Skip unstructured data for memory efficiency
            sparse_matrix=True,
            density=density,
            chunk_size=chunk_size,
            mode=mode,
            use_chunked_layers=False,  # Disable chunked layer processing
            minimal_metadata=True      # Use minimal metadata for large datasets
        )
    
    def create_extreme_dataset(self, n_cells: int = 1000000, n_genes: int = 5000,
                             density: float = 0.001, chunk_size: int = 50000,
                             mode: str = 'fast') -> sc.AnnData:
        """Create an extreme dataset (1M+ cells) for ultimate stress testing.
        
        This method is optimized for generating very large datasets efficiently.
        Use 'fast' or 'skeleton' mode for best performance.
        
        Args:
            n_cells: Number of cells (default: 1M)
            n_genes: Number of genes (default: 5K)
            density: Density of non-zero elements (default: 0.1%)
            chunk_size: Size of chunks for matrix generation
            mode: Generation mode - 'fast' or 'skeleton' recommended for large datasets
        """
        print(f"Creating extreme dataset: {n_cells:,} cells x {n_genes:,} genes")
        print(f"Mode: {mode}, Density: {density:.4f}")
        print(f"Estimated memory usage: {estimate_memory_usage(n_cells, n_genes, sparse=True):.1f}MB (sparse)")
        
        return self._create_anndata(
            n_cells=n_cells,
            n_genes=n_genes,
            add_missing=False,         # Skip missing values for speed
            add_dim_reductions=False,  # Skip dimensionality reductions
            add_layers=False,          # Skip layers for speed
            add_uns=False,             # Skip unstructured data
            sparse_matrix=True,
            density=density,
            chunk_size=chunk_size,
            mode=mode
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
                       density: float = 0.1,
                       chunk_size: int = 10000,
                       mode: str = 'realistic',
                       cell_types: Optional[List[str]] = None,
                       conditions: Optional[List[str]] = None,
                       gene_types: Optional[List[str]] = None,
                       use_chunked_layers: bool = False,
                       minimal_metadata: bool = False) -> sc.AnnData:
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
            n_cells, cell_types, conditions, add_missing, minimal_metadata
        )
        
        # Create gene metadata (var)
        var_data = self._create_var_data(
            n_genes, gene_types, add_missing, minimal_metadata
        )
        
        # Create expression matrix
        X = self._create_expression_matrix(n_cells, n_genes, sparse_matrix, density, chunk_size, mode)
        
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
                        conditions: List[str], add_missing: bool, minimal_metadata: bool = False) -> pd.DataFrame:
        """Create cell metadata DataFrame."""
        if minimal_metadata:
            # For large datasets, use minimal metadata to save memory
            obs_data = pd.DataFrame({
                'cell_type': pd.Categorical(
                    np.random.choice(cell_types, n_cells)
                ),
                'condition': pd.Categorical(
                    np.random.choice(conditions, n_cells)
                ),
                'quality_score': np.random.normal(0, 1, n_cells)
            })
        else:
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
                        add_missing: bool, minimal_metadata: bool = False) -> pd.DataFrame:
        """Create gene metadata DataFrame."""
        if minimal_metadata:
            # For large datasets, use minimal metadata to save memory
            var_data = pd.DataFrame({
                'gene_type': pd.Categorical(
                    np.random.choice(gene_types, n_genes)
                ),
                'chromosome': pd.Categorical(
                    [f'chr{i}' for i in np.random.randint(1, 23, n_genes)]
                ),
                'name': [f'GENE_{i:05d}' for i in range(n_genes)]
            })
        else:
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
                                sparse: bool = False, density: float = 0.1,
                                chunk_size: int = 10000, mode: str = 'realistic') -> Union[np.ndarray, scipy.sparse.spmatrix]:
        """Create expression matrix with realistic single-cell data patterns.
        
        Args:
            n_cells: Number of cells
            n_genes: Number of genes  
            sparse: Whether to create a sparse matrix
            density: Density of non-zero elements (0.0 to 1.0). Default 0.1 (10% non-zero)
            chunk_size: Size of chunks for large matrix generation
            mode: Generation mode - 'realistic' (unique indices), 'fast' (may have duplicates), 
                  or 'skeleton' (structure only, no values)
        """
        if sparse:
            # Ensure density is valid
            density = max(0.0, min(1.0, density))
            
            # For very large matrices, use chunked generation
            if n_cells * n_genes > 100_000_000:  # 100M elements threshold
                return self._create_chunked_sparse_matrix(n_cells, n_genes, density, chunk_size, mode)
            
            # For smaller matrices, use the existing optimized approach
            return self._create_single_sparse_matrix(n_cells, n_genes, density, mode)
        else:
            # Create dense matrix
            X = np.random.negative_binomial(5, 0.3, (n_cells, n_genes))
        
        return X
    
    def _create_single_sparse_matrix(self, n_cells: int, n_genes: int, density: float, mode: str) -> scipy.sparse.spmatrix:
        """Create a single sparse matrix using the optimized approach."""
        # Calculate number of non-zero elements
        nnz = int(n_cells * n_genes * density)
        
        # For very sparse matrices, use a more efficient approach
        if density < 0.01:  # Less than 1% density
            return self._create_very_sparse_matrix(n_cells, n_genes, nnz, mode)
        
        # For moderately sparse matrices, use optimized approach
        if density < 0.3:  # Less than 30% density
            return self._create_moderately_sparse_matrix(n_cells, n_genes, nnz, mode)
        
        # For dense matrices, use simple approach (duplicates are less likely)
        return self._create_dense_sparse_matrix(n_cells, n_genes, nnz, mode)
    
    def _create_chunked_sparse_matrix(self, n_cells: int, n_genes: int, density: float, 
                                    chunk_size: int, mode: str) -> scipy.sparse.spmatrix:
        """Create large sparse matrices using chunked generation."""
        import scipy.sparse as sp
        
        print(f"Generating chunked sparse matrix: {n_cells:,} cells × {n_genes:,} genes, density={density:.3f}")
        
        # Calculate total number of non-zero elements
        total_nnz = int(n_cells * n_genes * density)
        
        if mode == 'skeleton':
            # For skeleton mode, just create the structure without filling values
            print("Creating skeleton matrix (structure only)")
            return self._create_skeleton_matrix(n_cells, n_genes, total_nnz)
        
        # For very large matrices, use a more memory-efficient approach
        if total_nnz > 5_000_000:  # 5M non-zero elements threshold
            print("Using memory-efficient incremental CSR construction...")
            return self._create_incremental_csr_matrix(n_cells, n_genes, density, chunk_size, mode)
        
        # For smaller matrices, use the original approach
        all_rows = []
        all_cols = []
        all_values = []
        
        for chunk_start in range(0, n_cells, chunk_size):
            chunk_end = min(chunk_start + chunk_size, n_cells)
            chunk_cells = chunk_end - chunk_start
            
            # Adjust nnz for this chunk
            chunk_nnz = int(chunk_cells * n_genes * density)
            
            if mode == 'realistic':
                # Generate unique indices for this chunk
                chunk_rows, chunk_cols, chunk_vals = self._generate_unique_chunk_indices(
                    chunk_cells, n_genes, chunk_nnz, chunk_start
                )
            else:  # fast mode
                # Generate indices quickly (may have duplicates)
                chunk_rows, chunk_cols, chunk_vals = self._generate_fast_chunk_indices(
                    chunk_cells, n_genes, chunk_nnz, chunk_start
                )
            
            all_rows.extend(chunk_rows)
            all_cols.extend(chunk_cols)
            all_values.extend(chunk_vals)
            
            if chunk_start % (chunk_size * 10) == 0:
                print(f"  Processed {chunk_start:,}/{n_cells:,} cells")
        
        # Create the final sparse matrix
        print("Assembling final sparse matrix...")
        X = sp.csr_matrix((all_values, (all_rows, all_cols)), shape=(n_cells, n_genes))
        
        print(f"Final matrix: {X.shape}, nnz: {X.nnz:,}, density: {X.nnz/(n_cells*n_genes):.6f}")
        return X
    
    def _create_incremental_csr_matrix(self, n_cells: int, n_genes: int, density: float,
                                     chunk_size: int, mode: str) -> scipy.sparse.spmatrix:
        """Create large sparse matrices using incremental CSR construction to save memory."""
        import scipy.sparse as sp
        
        # Pre-allocate CSR matrix structure
        total_nnz = int(n_cells * n_genes * density)
        
        # Initialize CSR arrays
        indptr = np.zeros(n_cells + 1, dtype=np.int32)
        indices = np.zeros(total_nnz, dtype=np.int32)
        data = np.zeros(total_nnz, dtype=np.float32)
        
        current_nnz = 0
        
        for chunk_start in range(0, n_cells, chunk_size):
            chunk_end = min(chunk_start + chunk_size, n_cells)
            chunk_cells = chunk_end - chunk_start
            
            # Adjust nnz for this chunk
            chunk_nnz = int(chunk_cells * n_genes * density)
            
            if mode == 'realistic':
                # Generate unique indices for this chunk
                chunk_rows, chunk_cols, chunk_vals = self._generate_unique_chunk_indices(
                    chunk_cells, n_genes, chunk_nnz, chunk_start
                )
            else:  # fast mode
                # Generate indices quickly (may have duplicates)
                chunk_rows, chunk_cols, chunk_vals = self._generate_fast_chunk_indices(
                    chunk_cells, n_genes, chunk_nnz, chunk_start
                )
            
            # Sort by row for CSR format
            sorted_indices = np.argsort(chunk_rows)
            chunk_rows = chunk_rows[sorted_indices]
            chunk_cols = chunk_cols[sorted_indices]
            chunk_vals = chunk_vals[sorted_indices]
            
            # Fill CSR arrays
            for i, (row, col, val) in enumerate(zip(chunk_rows, chunk_cols, chunk_vals)):
                if current_nnz >= total_nnz:
                    break
                indices[current_nnz] = col
                data[current_nnz] = val
                current_nnz += 1
            
            # Update indptr for this chunk
            for row in range(chunk_start, chunk_end):
                indptr[row + 1] = current_nnz
            
            if chunk_start % (chunk_size * 10) == 0:
                print(f"  Processed {chunk_start:,}/{n_cells:,} cells")
        
        # Create CSR matrix
        X = sp.csr_matrix((data[:current_nnz], indices[:current_nnz], indptr), 
                         shape=(n_cells, n_genes))
        
        print(f"Final matrix: {X.shape}, nnz: {X.nnz:,}, density: {X.nnz/(n_cells*n_genes):.6f}")
        return X
    
    def _generate_unique_chunk_indices(self, chunk_cells: int, n_genes: int, nnz: int, 
                                     row_offset: int) -> tuple:
        """Generate unique indices for a chunk."""
        # Use rejection sampling for unique indices
        target_size = int(nnz * 1.2)  # Generate 20% more to account for duplicates
        
        # Generate initial indices
        rows = np.random.randint(0, chunk_cells, target_size) + row_offset
        cols = np.random.randint(0, n_genes, target_size)
        
        # Create unique pairs
        pairs = set(zip(rows, cols))
        
        # If we don't have enough unique pairs, generate more
        while len(pairs) < nnz:
            additional_size = min(1000, nnz - len(pairs))
            new_rows = np.random.randint(0, chunk_cells, additional_size) + row_offset
            new_cols = np.random.randint(0, n_genes, additional_size)
            pairs.update(zip(new_rows, new_cols))
        
        # Take exactly nnz pairs
        unique_pairs = list(pairs)[:nnz]
        rows = np.array([p[0] for p in unique_pairs])
        cols = np.array([p[1] for p in unique_pairs])
        values = self._generate_realistic_expression_values(nnz)
        
        return rows, cols, values
    
    def _generate_fast_chunk_indices(self, chunk_cells: int, n_genes: int, nnz: int, 
                                   row_offset: int) -> tuple:
        """Generate indices quickly (may have duplicates)."""
        rows = np.random.randint(0, chunk_cells, nnz) + row_offset
        cols = np.random.randint(0, n_genes, nnz)
        values = self._generate_realistic_expression_values(nnz)
        
        return rows, cols, values
    
    def _create_skeleton_matrix(self, n_cells: int, n_genes: int, nnz: int) -> scipy.sparse.spmatrix:
        """Create a skeleton matrix with structure but no meaningful values."""
        import scipy.sparse as sp
        
        # Generate random indices (no need to worry about duplicates for skeleton)
        rows = np.random.randint(0, n_cells, nnz)
        cols = np.random.randint(0, n_genes, nnz)
        
        # Use placeholder values (all 1s or random small integers)
        values = np.ones(nnz, dtype=np.int32)
        
        return sp.csr_matrix((values, (rows, cols)), shape=(n_cells, n_genes))
    
    def _create_moderately_sparse_matrix(self, n_cells: int, n_genes: int, nnz: int, mode: str) -> scipy.sparse.spmatrix:
        """Create moderately sparse matrices efficiently with guaranteed unique indices."""
        import scipy.sparse as sp
        
        if mode == 'fast':
            # Fast mode - may have duplicates
            cell_indices = np.random.randint(0, n_cells, nnz)
            gene_indices = np.random.randint(0, n_genes, nnz)
            values = self._generate_realistic_expression_values(nnz)
        elif mode == 'skeleton':
            # Skeleton mode - just structure
            cell_indices = np.random.randint(0, n_cells, nnz)
            gene_indices = np.random.randint(0, n_genes, nnz)
            values = np.ones(nnz, dtype=np.int32)
        else:  # realistic mode
            # Use a more efficient approach for moderate sparsity
            # Generate indices using rejection sampling with larger initial sample
            target_size = int(nnz * 1.2)  # Generate 20% more to account for duplicates
            
            # Generate initial indices
            cell_indices = np.random.randint(0, n_cells, target_size)
            gene_indices = np.random.randint(0, n_genes, target_size)
            
            # Create unique pairs
            pairs = set(zip(cell_indices, gene_indices))
            
            # If we don't have enough unique pairs, generate more
            while len(pairs) < nnz:
                additional_size = min(1000, nnz - len(pairs))
                new_cells = np.random.randint(0, n_cells, additional_size)
                new_genes = np.random.randint(0, n_genes, additional_size)
                pairs.update(zip(new_cells, new_genes))
            
            # Take exactly nnz pairs
            unique_pairs = list(pairs)[:nnz]
            cell_indices = np.array([p[0] for p in unique_pairs])
            gene_indices = np.array([p[1] for p in unique_pairs])
            values = self._generate_realistic_expression_values(nnz)
        
        # Create sparse matrix
        return sp.csr_matrix((values, (cell_indices, gene_indices)), 
                           shape=(n_cells, n_genes))
    
    def _create_dense_sparse_matrix(self, n_cells: int, n_genes: int, nnz: int, mode: str) -> scipy.sparse.spmatrix:
        """Create dense sparse matrices using simple approach (duplicates are less likely)."""
        import scipy.sparse as sp
        
        # For dense matrices, duplicates are less likely, so use simple approach
        cell_indices = np.random.randint(0, n_cells, nnz)
        gene_indices = np.random.randint(0, n_genes, nnz)
        
        if mode == 'skeleton':
            values = np.ones(nnz, dtype=np.int32)
        else:
            values = self._generate_realistic_expression_values(nnz)
        
        # Create sparse matrix
        return sp.csr_matrix((values, (cell_indices, gene_indices)), 
                           shape=(n_cells, n_genes))
    
    def _create_very_sparse_matrix(self, n_cells: int, n_genes: int, nnz: int, mode: str) -> Any:
        """Create very sparse matrices efficiently using COO format."""
        import scipy.sparse as sp
        
        # For very sparse matrices, use coordinate format for efficiency
        cell_indices = np.random.choice(n_cells, nnz, replace=True)
        gene_indices = np.random.choice(n_genes, nnz, replace=True)
        
        if mode == 'skeleton':
            values = np.ones(nnz, dtype=np.int32)
        else:
            values = self._generate_realistic_expression_values(nnz)
        
        # Create COO matrix and convert to CSR
        coo_matrix = sp.coo_matrix((values, (cell_indices, gene_indices)), 
                                  shape=(n_cells, n_genes))
        csr_matrix = coo_matrix.tocsr()
        return csr_matrix
    
    def _generate_realistic_expression_values(self, n_values: int) -> np.ndarray:
        """Generate realistic single-cell expression values.
        
        Single-cell data typically follows a negative binomial distribution
        with many zeros and a long tail of high expression values.
        """
        # Use negative binomial distribution for realistic counts
        # Parameters tuned for typical single-cell data
        r, p = 5, 0.3  # Shape and probability parameters
        
        # Generate base values
        values = np.random.negative_binomial(r, p, n_values)
        
        # Add some zeros to make it more realistic (some genes are truly not expressed)
        zero_prob = 0.3  # 30% chance of being zero
        zero_mask = np.random.random(n_values) < zero_prob
        values[zero_mask] = 0
        
        # Ensure no negative values (shouldn't happen with negative binomial, but just in case)
        values = np.maximum(values, 0)
        
        return values
    
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
                        add_missing: bool = True, density: float = 0.1) -> sc.AnnData:
    """Create large AnnData object for stress testing."""
    factory = MockAnnDataFactory()
    return factory.create_large(n_cells, n_genes, add_missing, density)


def create_memory_efficient_large_anndata(n_cells: int = 10000, n_genes: int = 5000,
                                         add_missing: bool = True, density: float = 0.1) -> sc.AnnData:
    """Create memory-efficient large AnnData object for stress testing.
    
    This function creates large datasets without dense layers to avoid
    excessive memory consumption during stress testing.
    """
    factory = MockAnnDataFactory()
    return factory.create_memory_efficient_large(n_cells, n_genes, add_missing, density)


def create_massive_anndata(n_cells: int = 100000, n_genes: int = 10000,
                          add_missing: bool = True, density: float = 0.1,
                          chunk_size: int = 10000, mode: str = 'realistic') -> sc.AnnData:
    """Create massive AnnData object for extreme stress testing.
    
    This function creates datasets with 100k+ cells using chunked operations
    to handle memory efficiently. Suitable for testing with real-world scale data.
    """
    factory = MockAnnDataFactory()
    return factory.create_massive_dataset(n_cells, n_genes, add_missing, density, chunk_size, mode)


def create_extreme_anndata(n_cells: int = 1000000, n_genes: int = 5000,
                          density: float = 0.001, chunk_size: int = 50000,
                          mode: str = 'fast') -> sc.AnnData:
    """Create extreme AnnData object for ultimate stress testing.
    
    This function creates datasets with 1M+ cells using optimized chunked operations.
    Use 'fast' or 'skeleton' mode for best performance with large datasets.
    """
    factory = MockAnnDataFactory()
    return factory.create_extreme_dataset(n_cells, n_genes, density, chunk_size, mode)


def create_fast_large_anndata(n_cells: int = 100000, n_genes: int = 5000,
                             density: float = 0.1, chunk_size: int = 10000) -> sc.AnnData:
    """Create large AnnData object using fast generation mode.
    
    This function prioritizes speed over perfect accuracy (may have duplicate indices).
    Suitable for stress testing where speed is more important than data quality.
    """
    factory = MockAnnDataFactory()
    return factory.create_massive_dataset(n_cells, n_genes, add_missing=False, 
                                        density=density, chunk_size=chunk_size, mode='fast')


def create_skeleton_anndata(n_cells: int = 100000, n_genes: int = 5000,
                           density: float = 0.1, chunk_size: int = 10000) -> sc.AnnData:
    """Create large AnnData object with skeleton matrix (structure only).
    
    This function creates a matrix with the correct structure but placeholder values.
    Fastest option for testing pipeline structure without realistic data.
    """
    factory = MockAnnDataFactory()
    return factory.create_massive_dataset(n_cells, n_genes, add_missing=False,
                                        density=density, chunk_size=chunk_size, mode='skeleton')


def create_edge_case_anndata() -> sc.AnnData:
    """Create AnnData object with various edge cases and problematic data."""
    factory = MockAnnDataFactory()
    return factory.create_edge_cases()


# Utility functions for testing
def _check_missing_values(data):
    """Check for missing values in a Dataset2D or DataFrame object."""
    try:
        if hasattr(data, 'to_pandas'):
            df = data.to_pandas()  # type: ignore
        elif hasattr(data, 'to_dataframe'):
            df = data.to_dataframe()  # type: ignore
        else:
            # Assume it's already a DataFrame or similar
            df = data
        return bool(df.isnull().values.any())
    except (AttributeError, TypeError):
        # If we can't check for missing values, assume False
        return False

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
        'has_missing_obs': _check_missing_values(adata.obs),
        'has_missing_var': _check_missing_values(adata.var),
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