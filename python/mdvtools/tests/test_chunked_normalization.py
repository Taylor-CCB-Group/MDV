#!/usr/bin/env python3
"""
Test chunked normalization functionality for large datasets.

This module tests the chunked normalization utilities that allow processing
of massive single-cell datasets without memory explosions.
"""

import numpy as np
import scipy.sparse
import pytest
import time
import psutil
import gc
from typing import Dict

from .mock_anndata import (
    chunked_log1p_normalization,
    chunked_zscore_normalization,
    estimate_memory_usage
)


def get_memory_usage() -> Dict[str, float]:
    """Get current memory usage in MB."""
    process = psutil.Process()
    memory_info = process.memory_info()
    return {
        'rss': memory_info.rss / 1024 / 1024,  # Resident Set Size in MB
        'vms': memory_info.vms / 1024 / 1024,  # Virtual Memory Size in MB
    }


class TestChunkedNormalization:
    """Test class for chunked normalization utilities."""
    
    def test_chunked_log1p_small_matrix(self):
        """Test chunked log1p normalization on small matrices."""
        # Create small sparse matrix
        n_cells, n_genes = 100, 50
        sparsity = 0.9
        nnz = int(n_cells * n_genes * (1 - sparsity))
        
        # Generate random data
        cell_indices = np.random.randint(0, n_cells, nnz)
        gene_indices = np.random.randint(0, n_genes, nnz)
        values = np.random.negative_binomial(5, 0.3, nnz)
        
        sparse_matrix = scipy.sparse.csr_matrix(
            (values, (cell_indices, gene_indices)), 
            shape=(n_cells, n_genes)
        )
        
        # Test chunked normalization
        result = chunked_log1p_normalization(sparse_matrix, chunk_size=50)
        
        # Verify result
        assert scipy.sparse.issparse(result)
        assert result.shape == (n_cells, n_genes)
        
        # Compare with direct computation
        expected = scipy.sparse.csc_matrix(np.log1p(sparse_matrix.toarray()))
        np.testing.assert_array_almost_equal(result.toarray(), expected.toarray())
    
    def test_chunked_log1p_large_matrix(self):
        """Test chunked log1p normalization on large matrices."""
        # Create large sparse matrix
        n_cells, n_genes = 10000, 5000
        sparsity = 0.95
        nnz = int(n_cells * n_genes * (1 - sparsity))
        
        print(f"Creating large sparse matrix: {n_cells:,} cells x {n_genes:,} genes")
        print(f"Estimated memory: {estimate_memory_usage(n_cells, n_genes, sparse=True):.1f}MB")
        
        # Generate random data
        cell_indices = np.random.randint(0, n_cells, nnz)
        gene_indices = np.random.randint(0, n_genes, nnz)
        values = np.random.negative_binomial(5, 0.3, nnz)
        
        sparse_matrix = scipy.sparse.csr_matrix(
            (values, (cell_indices, gene_indices)), 
            shape=(n_cells, n_genes)
        )
        
        # Test chunked normalization
        start_time = time.time()
        start_memory = get_memory_usage()
        
        result = chunked_log1p_normalization(sparse_matrix, chunk_size=1000)
        
        end_time = time.time()
        end_memory = get_memory_usage()
        
        print(f"Chunked log1p time: {end_time - start_time:.2f}s")
        print(f"Memory increase: {end_memory['rss'] - start_memory['rss']:.1f}MB")
        
        # Verify result
        assert scipy.sparse.issparse(result)
        assert result.shape == (n_cells, n_genes)
        
        # Check that we didn't use excessive memory
        memory_increase = end_memory['rss'] - start_memory['rss']
        assert memory_increase < 2000, f"Memory increase too high: {memory_increase:.1f}MB"
    
    def test_chunked_zscore_small_matrix(self):
        """Test chunked z-score normalization on small matrices."""
        # Create small sparse matrix
        n_cells, n_genes = 100, 50
        sparsity = 0.9
        nnz = int(n_cells * n_genes * (1 - sparsity))
        
        # Generate random data
        cell_indices = np.random.randint(0, n_cells, nnz)
        gene_indices = np.random.randint(0, n_genes, nnz)
        values = np.random.negative_binomial(5, 0.3, nnz)
        
        sparse_matrix = scipy.sparse.csr_matrix(
            (values, (cell_indices, gene_indices)), 
            shape=(n_cells, n_genes)
        )
        
        # Test chunked normalization
        result = chunked_zscore_normalization(sparse_matrix, chunk_size=50)
        
        # Verify result
        assert scipy.sparse.issparse(result)
        assert result.shape == (n_cells, n_genes)
        
        # Compare with direct computation
        dense_matrix = sparse_matrix.toarray()
        expected = (dense_matrix - dense_matrix.mean(axis=0)) / (dense_matrix.std(axis=0) + 1e-8)
        np.testing.assert_array_almost_equal(result.toarray(), expected, decimal=5)
    
    def test_chunked_zscore_large_matrix(self):
        """Test chunked z-score normalization on large matrices."""
        # Create large sparse matrix
        n_cells, n_genes = 10000, 5000
        sparsity = 0.95
        nnz = int(n_cells * n_genes * (1 - sparsity))
        
        print(f"Creating large sparse matrix for z-score: {n_cells:,} cells x {n_genes:,} genes")
        
        # Generate random data
        cell_indices = np.random.randint(0, n_cells, nnz)
        gene_indices = np.random.randint(0, n_genes, nnz)
        values = np.random.negative_binomial(5, 0.3, nnz)
        
        sparse_matrix = scipy.sparse.csr_matrix(
            (values, (cell_indices, gene_indices)), 
            shape=(n_cells, n_genes)
        )
        
        # Test chunked normalization
        start_time = time.time()
        start_memory = get_memory_usage()
        
        result = chunked_zscore_normalization(sparse_matrix, chunk_size=1000)
        
        end_time = time.time()
        end_memory = get_memory_usage()
        
        print(f"Chunked z-score time: {end_time - start_time:.2f}s")
        print(f"Memory increase: {end_memory['rss'] - start_memory['rss']:.1f}MB")
        
        # Verify result
        assert scipy.sparse.issparse(result)
        assert result.shape == (n_cells, n_genes)
        
        # Check that we didn't use excessive memory
        memory_increase = end_memory['rss'] - start_memory['rss']
        assert memory_increase < 3000, f"Memory increase too high: {memory_increase:.1f}MB"
    
    def test_memory_estimation(self):
        """Test memory usage estimation function."""
        # Test sparse matrix estimation
        sparse_memory = estimate_memory_usage(100000, 10000, sparse=True)
        print(f"Sparse 100k x 10k matrix: {sparse_memory:.1f}MB")
        assert sparse_memory > 0
        assert sparse_memory < 10000  # Should be reasonable
        
        # Test dense matrix estimation
        dense_memory = estimate_memory_usage(100000, 10000, sparse=False)
        print(f"Dense 100k x 10k matrix: {dense_memory:.1f}MB")
        assert dense_memory > 0
        assert dense_memory > sparse_memory  # Dense should use more memory
    
    @pytest.mark.performance
    def test_chunked_vs_traditional_performance(self):
        """Compare performance of chunked vs traditional normalization."""
        # Create moderately large matrix
        n_cells, n_genes = 50000, 5000
        sparsity = 0.95
        nnz = int(n_cells * n_genes * (1 - sparsity))
        
        print("\nComparing chunked vs traditional normalization")
        print(f"Matrix size: {n_cells:,} cells x {n_genes:,} genes")
        
        # Generate random data
        cell_indices = np.random.randint(0, n_cells, nnz)
        gene_indices = np.random.randint(0, n_genes, nnz)
        values = np.random.negative_binomial(5, 0.3, nnz)
        
        sparse_matrix = scipy.sparse.csr_matrix(
            (values, (cell_indices, gene_indices)), 
            shape=(n_cells, n_genes)
        )
        
        # Test chunked approach
        gc.collect()
        start_time = time.time()
        start_memory = get_memory_usage()
        
        chunked_result = chunked_log1p_normalization(sparse_matrix, chunk_size=1000)
        
        chunked_time = time.time() - start_time
        chunked_memory = get_memory_usage()
        chunked_memory_increase = chunked_memory['rss'] - start_memory['rss']
        
        print(f"Chunked approach: {chunked_time:.2f}s, Memory: {chunked_memory_increase:.1f}MB")
        
        # Verify chunked result is reasonable
        assert scipy.sparse.issparse(chunked_result)
        assert chunked_result.shape == (n_cells, n_genes)
        assert chunked_memory_increase < 2000  # Should be memory efficient


if __name__ == "__main__":
    pytest.main([__file__]) 