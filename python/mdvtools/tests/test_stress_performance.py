#!/usr/bin/env python3
"""
Stress testing and performance testing for MDV conversion pipeline.

This module contains tests that verify the conversion pipeline can handle
large datasets, edge cases, and various data configurations efficiently.
"""

import os
import tempfile
import shutil
import time
import pytest
import psutil
import gc
import pandas as pd
import scipy.sparse
from typing import Dict
from contextlib import contextmanager

from mdvtools.conversions import convert_scanpy_to_mdv
from mdvtools.mdvproject import MDVProject
from .mock_anndata import (
    MockAnnDataFactory,
    suppress_anndata_warnings,
)


def get_memory_usage() -> Dict[str, float]:
    """Get current memory usage in MB."""
    process = psutil.Process(os.getpid())
    memory_info = process.memory_info()
    return {
        'rss': memory_info.rss / 1024 / 1024,  # Resident Set Size in MB
        'vms': memory_info.vms / 1024 / 1024,  # Virtual Memory Size in MB
    }


@contextmanager
def temp_mdv_project():
    """Context manager for temporary MDV project creation and cleanup."""
    test_dir = tempfile.mkdtemp()
    try:
        yield test_dir
    finally:
        if os.path.exists(test_dir):
            shutil.rmtree(test_dir)


class TestStressTesting:
    """Test class for stress testing the MDV conversion pipeline."""
    
    @pytest.mark.performance
    def test_large_dataset_conversion(self):
        """Test conversion of large datasets (20k cells, 3k genes)."""
        factory = MockAnnDataFactory(random_seed=42)
        
        # Create large dataset
        start_time = time.time()
        # start_memory = get_memory_usage()
        
        adata = factory.create_large(20000, 3000)
        
        creation_time = time.time() - start_time
        creation_memory = get_memory_usage()
        
        print(f"Large dataset creation: {creation_time:.2f}s, Memory: {creation_memory['rss']:.1f}MB")
        
        # Test conversion
        with temp_mdv_project() as test_dir:
            with suppress_anndata_warnings():
                conversion_start = time.time()
                mdv = convert_scanpy_to_mdv(test_dir, adata, delete_existing=True, chunk_data=True)
                conversion_time = time.time() - conversion_start
                
                conversion_memory = get_memory_usage()
                print(f"Conversion time: {conversion_time:.2f}s, Memory: {conversion_memory['rss']:.1f}MB")
        
        # Verify results
        assert isinstance(mdv, MDVProject)
        assert "cells" in mdv.get_datasource_names()
        assert "genes" in mdv.get_datasource_names()
        
        # Check performance metrics
        assert conversion_time < 120, f"Conversion took too long: {conversion_time:.2f}s"
        assert conversion_memory['rss'] < 4000, f"Memory usage too high: {conversion_memory['rss']:.1f}MB"
    
    @pytest.mark.performance
    def test_very_large_dataset_conversion(self):
        """Test conversion of very large datasets (50k cells, 5k genes)."""
        factory = MockAnnDataFactory(random_seed=42)
        
        # Create very large dataset with memory-efficient method
        start_time = time.time()
        adata = factory.create_memory_efficient_large(50000, 5000)
        creation_time = time.time() - start_time
        
        print(f"Very large dataset creation: {creation_time:.2f}s")
        print(f"Dataset size: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
        
        # Test conversion with chunking
        with temp_mdv_project() as test_dir:
            with suppress_anndata_warnings():
                conversion_start = time.time()
                mdv = convert_scanpy_to_mdv(
                    test_dir, adata, delete_existing=True, chunk_data=True
                )
                conversion_time = time.time() - conversion_start
                
                print(f"Very large conversion time: {conversion_time:.2f}s")
        
        # Verify results
        assert isinstance(mdv, MDVProject)
        assert "cells" in mdv.get_datasource_names()
        assert "genes" in mdv.get_datasource_names()
        
        # Check performance metrics
        assert conversion_time < 180, f"Conversion took too long: {conversion_time:.2f}s"
    
    @pytest.mark.performance
    def test_massive_dataset_conversion(self):
        """Test conversion of massive datasets (100k cells, 8k genes)."""
        factory = MockAnnDataFactory(random_seed=42)
        
        # Create massive dataset using chunked operations
        start_time = time.time()
        adata = factory.create_massive_dataset(100000, 8000)
        creation_time = time.time() - start_time
        
        print(f"Massive dataset creation: {creation_time:.2f}s")
        print(f"Dataset size: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
        
        # Test conversion with chunking
        with temp_mdv_project() as test_dir:
            with suppress_anndata_warnings():
                conversion_start = time.time()
                mdv = convert_scanpy_to_mdv(
                    test_dir, adata, delete_existing=True, chunk_data=True
                )
                conversion_time = time.time() - conversion_start
                
                print(f"Massive conversion time: {conversion_time:.2f}s")
        
        # Verify results
        assert isinstance(mdv, MDVProject)
        assert "cells" in mdv.get_datasource_names()
        assert "genes" in mdv.get_datasource_names()
        
        # Check performance metrics
        assert conversion_time < 300, f"Conversion took too long: {conversion_time:.2f}s"
    
    @pytest.mark.performance
    def test_memory_efficiency(self):
        """Test memory efficiency during conversion."""
        factory = MockAnnDataFactory(random_seed=42)
        
        # Create dataset with reasonable size
        adata = factory.create_realistic(10000, 5000)
        
        # Force garbage collection before conversion
        gc.collect()
        initial_memory = get_memory_usage()
        
        with temp_mdv_project() as test_dir:
            with suppress_anndata_warnings():
                mdv = convert_scanpy_to_mdv(test_dir, adata, delete_existing=True)
                
                # Force garbage collection after conversion
                gc.collect()
                final_memory = get_memory_usage()
                
                memory_increase = final_memory['rss'] - initial_memory['rss']
                print(f"Memory increase during conversion: {memory_increase:.1f}MB")
        
        # Verify conversion worked
        assert isinstance(mdv, MDVProject)
        
        # Check memory efficiency (should not increase by more than 1GB)
        assert memory_increase < 1000, f"Memory increase too high: {memory_increase:.1f}MB"
    
    @pytest.mark.performance
    def test_multiple_conversions_same_session(self):
        """Test multiple conversions in the same session to check for memory leaks."""
        factory = MockAnnDataFactory(random_seed=42)
        
        # Force garbage collection
        gc.collect()
        initial_memory = get_memory_usage()
        
        for i in range(3):  # Reduced from 5 to 3 for larger datasets
            # Create new dataset each time with smaller size
            adata = factory.create_realistic(5000, 2000)
            
            with temp_mdv_project() as test_dir:
                with suppress_anndata_warnings():
                    convert_scanpy_to_mdv(test_dir, adata, delete_existing=True)
            
            # Force garbage collection
            gc.collect()
            current_memory = get_memory_usage()
            
            print(f"Iteration {i+1}: Memory usage: {current_memory['rss']:.1f}MB")
        
        final_memory = get_memory_usage()
        memory_increase = final_memory['rss'] - initial_memory['rss']
        
        # Check for memory leaks (should not increase by more than 500MB over 3 iterations)
        assert memory_increase < 500, f"Potential memory leak: {memory_increase:.1f}MB increase over 3 iterations"
    
    @pytest.mark.performance
    def test_chunked_normalization_performance(self):
        """Test performance of chunked normalization vs traditional approach."""
        factory = MockAnnDataFactory(random_seed=42)
        
        # Test with a moderately large dataset
        n_cells, n_genes = 25000, 3000
        
        print(f"\nTesting chunked normalization with {n_cells:,} cells x {n_genes:,} genes")
        
        # Test memory-efficient approach
        start_time = time.time()
        # start_memory = get_memory_usage()
        
        adata_chunked = factory.create_massive_dataset(n_cells, n_genes)
        
        creation_time = time.time() - start_time
        creation_memory = get_memory_usage()
        
        print(f"Chunked creation: {creation_time:.2f}s, Memory: {creation_memory['rss']:.1f}MB")
        
        # Verify chunked layers were created
        assert 'log1p' in adata_chunked.layers
        assert 'scaled' in adata_chunked.layers
        
        # Test conversion
        with temp_mdv_project() as test_dir:
            with suppress_anndata_warnings():
                conversion_start = time.time()
                mdv = convert_scanpy_to_mdv(test_dir, adata_chunked, delete_existing=True, chunk_data=True)
                conversion_time = time.time() - conversion_start
        
        print(f"Chunked conversion: {conversion_time:.2f}s")
        
        # Verify results
        assert isinstance(mdv, MDVProject)
        assert "cells" in mdv.get_datasource_names()
        assert "genes" in mdv.get_datasource_names()


class TestEdgeCaseStressTesting:
    """Test class for edge case stress testing."""
    
    @pytest.mark.performance
    def test_edge_case_conversion(self):
        """Test conversion of edge case data."""
        factory = MockAnnDataFactory(random_seed=42)
        
        adata = factory.create_edge_cases()
        
        with temp_mdv_project() as test_dir:
            with suppress_anndata_warnings():
                mdv = convert_scanpy_to_mdv(test_dir, adata, delete_existing=True)
        
        # Verify conversion worked despite edge cases
        assert isinstance(mdv, MDVProject)
        assert "cells" in mdv.get_datasource_names()
        assert "genes" in mdv.get_datasource_names()
    
    @pytest.mark.performance
    def test_mixed_data_types(self):
        """Test conversion with mixed data types."""
        factory = MockAnnDataFactory(random_seed=42)
        
        # Create dataset with various data types
        adata = factory.create_with_specific_features(
            n_cells=100,
            n_genes=200,
            add_missing=True
        )
        
        # Add some additional mixed data types
        adata.obs['mixed_column'] = [1, 'string', True, 3.14, None] * 20
        adata.var['mixed_column'] = [1, 'string', True, 3.14, None] * 40
        
        with temp_mdv_project() as test_dir:
            with suppress_anndata_warnings():
                mdv = convert_scanpy_to_mdv(test_dir, adata, delete_existing=True)
        
        # Verify conversion worked
        assert isinstance(mdv, MDVProject)
        assert "cells" in mdv.get_datasource_names()
        assert "genes" in mdv.get_datasource_names()
    
    @pytest.mark.performance
    def test_extreme_categorical_values(self):
        """Test conversion with extreme categorical values."""
        factory = MockAnnDataFactory(random_seed=42)
        
        # Create dataset with extreme categorical values
        adata = factory.create_minimal(100, 50)
        
        # Add extreme categorical values
        adata.obs['extreme_cat'] = pd.Categorical([''] * 50 + ['A' * 1000] * 50)
        adata.var['extreme_cat'] = pd.Categorical([''] * 25 + ['B' * 1000] * 25)
        
        with temp_mdv_project() as test_dir:
            with suppress_anndata_warnings():
                mdv = convert_scanpy_to_mdv(test_dir, adata, delete_existing=True)
        
        # Verify conversion worked
        assert isinstance(mdv, MDVProject)
        assert "cells" in mdv.get_datasource_names()
        assert "genes" in mdv.get_datasource_names()


class TestPerformanceBenchmarks:
    """Test class for performance benchmarking."""
    
    @pytest.mark.performance
    def test_conversion_speed_benchmark(self):
        """Benchmark conversion speed for different dataset sizes."""
        factory = MockAnnDataFactory(random_seed=42)
        
        sizes = [
            (1000, 500),
            (5000, 2000),
            (10000, 5000),
            (25000, 8000),
            (50000, 12000)
        ]
        
        results = {}
        
        for n_cells, n_genes in sizes:
            print(f"\nTesting size: {n_cells:,} cells x {n_genes:,} genes")
            
            # Create dataset
            start_time = time.time()
            adata = factory.create_massive_dataset(n_cells, n_genes)
            creation_time = time.time() - start_time
            
            # Convert dataset
            with temp_mdv_project() as test_dir:
                with suppress_anndata_warnings():
                    conversion_start = time.time()
                    convert_scanpy_to_mdv(test_dir, adata, delete_existing=True, chunk_data=True)
                    conversion_time = time.time() - conversion_start
            
            # Calculate metrics
            total_time = creation_time + conversion_time
            cells_per_second = n_cells / conversion_time if conversion_time > 0 else 0
            genes_per_second = n_genes / conversion_time if conversion_time > 0 else 0
            
            results[(n_cells, n_genes)] = {
                'creation_time': creation_time,
                'conversion_time': conversion_time,
                'total_time': total_time,
                'cells_per_second': cells_per_second,
                'genes_per_second': genes_per_second
            }
            
            print(f"  Creation: {creation_time:.2f}s")
            print(f"  Conversion: {conversion_time:.2f}s")
            print(f"  Total: {total_time:.2f}s")
            print(f"  Rate: {cells_per_second:.0f} cells/s, {genes_per_second:.0f} genes/s")
        
        # Verify all conversions worked
        for (n_cells, n_genes), metrics in results.items():
            assert metrics['conversion_time'] > 0, f"Conversion failed for {n_cells:,}x{n_genes:,}"
            assert metrics['cells_per_second'] > 0, f"Invalid conversion rate for {n_cells:,}x{n_genes:,}"
    
    @pytest.mark.performance
    def test_memory_usage_benchmark(self):
        """Benchmark memory usage for different dataset sizes."""
        factory = MockAnnDataFactory(random_seed=42)
        
        sizes = [
            (1000, 500),
            (5000, 2000),
            (10000, 5000)
        ]
        
        results = {}
        
        for n_cells, n_genes in sizes:
            print(f"\nTesting memory for: {n_cells:,} cells x {n_genes:,} genes")
            
            # Force garbage collection
            gc.collect()
            initial_memory = get_memory_usage()
            
            # Create and convert dataset
            adata = factory.create_massive_dataset(n_cells, n_genes)
            
            with temp_mdv_project() as test_dir:
                with suppress_anndata_warnings():
                    convert_scanpy_to_mdv(test_dir, adata, delete_existing=True, chunk_data=True)
            
            # Force garbage collection
            gc.collect()
            final_memory = get_memory_usage()
            
            memory_increase = final_memory['rss'] - initial_memory['rss']
            memory_per_cell = memory_increase / n_cells if n_cells > 0 else 0
            memory_per_gene = memory_increase / n_genes if n_genes > 0 else 0
            
            results[(n_cells, n_genes)] = {
                'memory_increase': memory_increase,
                'memory_per_cell': memory_per_cell,
                'memory_per_gene': memory_per_gene
            }
            
            print(f"  Memory increase: {memory_increase:.1f}MB")
            print(f"  Per cell: {memory_per_cell:.3f}MB")
            print(f"  Per gene: {memory_per_gene:.3f}MB")
        
        # Verify reasonable memory usage
        for (n_cells, n_genes), metrics in results.items():
            assert metrics['memory_increase'] >= 0, f"Negative memory increase for {n_cells:,}x{n_genes:,}"
            assert metrics['memory_per_cell'] < 0.1, f"Too much memory per cell for {n_cells:,}x{n_genes:,}"

    @pytest.mark.performance
    def test_memory_improvements_verification(self):
        """Test that memory improvements work without MDV project structure issues."""
        factory = MockAnnDataFactory(random_seed=42)
        
        # Test with a moderately large dataset
        n_cells, n_genes = 25000, 3000
        
        print(f"\nTesting memory improvements with {n_cells:,} cells x {n_genes:,} genes")
        
        # Test memory-efficient approach
        start_time = time.time()
        start_memory = get_memory_usage()
        
        adata_chunked = factory.create_massive_dataset(n_cells, n_genes)
        
        creation_time = time.time() - start_time
        creation_memory = get_memory_usage()
        
        print(f"Chunked creation: {creation_time:.2f}s, Memory: {creation_memory['rss']:.1f}MB")
        
        # Verify chunked layers were created
        assert 'log1p' in adata_chunked.layers
        assert 'scaled' in adata_chunked.layers
        
        # Verify the layers are sparse matrices (memory efficient)
        assert scipy.sparse.issparse(adata_chunked.layers['log1p'])
        assert scipy.sparse.issparse(adata_chunked.layers['scaled'])
        
        # Check memory usage is reasonable
        memory_increase = creation_memory['rss'] - start_memory['rss']
        print(f"Memory increase: {memory_increase:.1f}MB")
        
        # Should not use more than 2GB for this dataset size
        assert memory_increase < 2000, f"Memory usage too high: {memory_increase:.1f}MB"
        
        # Verify dataset properties
        assert adata_chunked.n_obs == n_cells
        assert adata_chunked.n_vars == n_genes
        assert scipy.sparse.issparse(adata_chunked.X)
        
        print("✓ Memory improvements verified successfully!")


if __name__ == "__main__":
    pytest.main([__file__]) 