#!/usr/bin/env python3
"""
Example usage of the mock AnnData module.

This script demonstrates how to use the MockAnnDataFactory and related utilities
for creating test data and stress testing the MDV conversion pipeline.
"""

import os
import tempfile
import shutil
from contextlib import contextmanager

from mdvtools.conversions import convert_scanpy_to_mdv
from mdvtools.mdvproject import MDVProject
from .mock_anndata import (
    MockAnnDataFactory,
    create_minimal_anndata,
    create_realistic_anndata,
    create_large_anndata,
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


def example_basic_usage():
    """Example of basic usage with convenience functions."""
    print("=== Basic Usage Examples ===")
    
    # Create minimal AnnData for quick testing
    print("\n1. Creating minimal AnnData...")
    adata_minimal = create_minimal_anndata(n_cells=20, n_genes=10)
    print(f"   Created: {adata_minimal.n_obs} cells x {adata_minimal.n_vars} genes")
    
    # Create realistic AnnData with typical features
    print("\n2. Creating realistic AnnData...")
    adata_realistic = create_realistic_anndata(n_cells=500, n_genes=1000)
    print(f"   Created: {adata_realistic.n_obs} cells x {adata_realistic.n_vars} genes")
    print(f"   Has PCA: {'X_pca' in adata_realistic.obsm}")
    print(f"   Has UMAP: {'X_umap' in adata_realistic.obsm}")
    print(f"   Has layers: {list(adata_realistic.layers.keys())}")
    
    # Create large AnnData for stress testing
    print("\n3. Creating large AnnData...")
    adata_large = create_large_anndata(n_cells=2000, n_genes=3000)
    print(f"   Created: {adata_large.n_obs} cells x {adata_large.n_vars} genes")
    print(f"   Is sparse: {hasattr(adata_large.X, 'toarray')}")
    
    # Create edge case AnnData
    print("\n4. Creating edge case AnnData...")
    adata_edge = create_edge_case_anndata()
    print(f"   Created: {adata_edge.n_obs} cells x {adata_edge.n_vars} genes")
    print(f"   Edge case columns: {list(adata_edge.obs.columns)}")


def example_factory_usage():
    """Example of using the MockAnnDataFactory class."""
    print("\n=== Factory Usage Examples ===")
    
    # Create factory with fixed random seed for reproducible results
    factory = MockAnnDataFactory(random_seed=42)
    
    # Create AnnData with specific features
    print("\n1. Creating AnnData with custom cell types...")
    custom_cell_types = ['Neuron', 'Astrocyte', 'Oligodendrocyte', 'Microglia']
    custom_conditions = ['Control', 'Disease', 'Treatment']
    
    adata_custom = factory.create_with_specific_features(
        cell_types=custom_cell_types,
        conditions=custom_conditions,
        n_cells=100,
        n_genes=200,
        add_missing=True
    )
    
    print(f"   Cell types: {list(adata_custom.obs['cell_type'].cat.categories)}")
    print(f"   Conditions: {list(adata_custom.obs['condition'].cat.categories)}")
    
    # Create different sizes for testing
    print("\n2. Creating datasets of different sizes...")
    sizes = [(50, 25), (200, 100), (1000, 500)]
    
    for n_cells, n_genes in sizes:
        adata = factory.create_realistic(n_cells, n_genes)
        summary = get_anndata_summary(adata)
        print(f"   {n_cells}x{n_genes}: {summary['n_cells']} cells, {summary['n_genes']} genes")


def example_conversion_testing():
    """Example of testing MDV conversion with mock data."""
    print("\n=== Conversion Testing Examples ===")
    
    factory = MockAnnDataFactory(random_seed=42)
    
    # Test conversion with minimal data
    print("\n1. Testing minimal data conversion...")
    adata_minimal = factory.create_minimal(50, 25)
    
    with temp_mdv_project() as test_dir:
        with suppress_anndata_warnings():
            mdv = convert_scanpy_to_mdv(test_dir, adata_minimal, delete_existing=True)
        
        print(f"   Conversion successful: {isinstance(mdv, MDVProject)}")
        print(f"   Datasources: {mdv.get_datasource_names()}")
    
    # Test conversion with realistic data
    print("\n2. Testing realistic data conversion...")
    adata_realistic = factory.create_realistic(200, 100)
    
    with temp_mdv_project() as test_dir:
        with suppress_anndata_warnings():
            mdv = convert_scanpy_to_mdv(test_dir, adata_realistic, delete_existing=True)
        
        # Check that all metadata was converted
        cells_metadata = mdv.get_datasource_metadata("cells")
        genes_metadata = mdv.get_datasource_metadata("genes")
        
        print(f"   Cells columns: {len(cells_metadata['columns'])}")
        print(f"   Genes columns: {len(genes_metadata['columns'])}")
    
    # Test conversion with edge cases
    print("\n3. Testing edge case conversion...")
    adata_edge = factory.create_edge_cases()
    
    with temp_mdv_project() as test_dir:
        with suppress_anndata_warnings():
            mdv = convert_scanpy_to_mdv(test_dir, adata_edge, delete_existing=True)
        
        print(f"   Edge case conversion successful: {isinstance(mdv, MDVProject)}")


def example_validation_and_summary():
    """Example of using validation and summary utilities."""
    print("\n=== Validation and Summary Examples ===")
    
    factory = MockAnnDataFactory(random_seed=42)
    
    # Create and validate different types of AnnData
    datasets = [
        ("minimal", factory.create_minimal(20, 10)),
        ("realistic", factory.create_realistic(100, 50)),
        ("large", factory.create_large(500, 200)),
        ("edge_cases", factory.create_edge_cases())
    ]
    
    for name, adata in datasets:
        print(f"\n{name.capitalize()} dataset:")
        
        # Validate structure
        is_valid = validate_anndata(adata)
        print(f"   Valid: {is_valid}")
        
        # Get summary
        summary = get_anndata_summary(adata)
        print(f"   Size: {summary['n_cells']} cells x {summary['n_genes']} genes")
        print(f"   Categorical obs: {len(summary['categorical_obs'])}")
        print(f"   Categorical var: {len(summary['categorical_var'])}")
        print(f"   Has missing obs: {summary['has_missing_obs']}")
        print(f"   Has missing var: {summary['has_missing_var']}")
        print(f"   Is sparse: {summary['sparse']}")


def example_stress_testing():
    """Example of stress testing with large datasets."""
    print("\n=== Stress Testing Examples ===")
    
    factory = MockAnnDataFactory(random_seed=42)
    
    # Test with progressively larger datasets
    test_sizes = [
        (100, 50, "small"),
        (500, 200, "medium"),
        (1000, 500, "large"),
        (2000, 1000, "very large")
    ]
    
    for n_cells, n_genes, size_name in test_sizes:
        print(f"\nTesting {size_name} dataset ({n_cells} cells x {n_genes} genes)...")
        
        # Create dataset
        adata = factory.create_realistic(n_cells, n_genes)
        
        # Test conversion
        with temp_mdv_project() as test_dir:
            with suppress_anndata_warnings():
                mdv = convert_scanpy_to_mdv(test_dir, adata, delete_existing=True)
            
            # Verify conversion
            cells_metadata = mdv.get_datasource_metadata("cells")
            genes_metadata = mdv.get_datasource_metadata("genes")
            
            print(f"   Conversion successful")
            print(f"   Cells columns: {len(cells_metadata['columns'])}")
            print(f"   Genes columns: {len(genes_metadata['columns'])}")


def main():
    """Run all examples."""
    print("Mock AnnData Module Usage Examples")
    print("=" * 50)
    
    try:
        example_basic_usage()
        example_factory_usage()
        example_conversion_testing()
        example_validation_and_summary()
        example_stress_testing()
        
        print("\n" + "=" * 50)
        print("All examples completed successfully!")
        
    except Exception as e:
        print(f"\nError running examples: {e}")
        raise


if __name__ == "__main__":
    main() 