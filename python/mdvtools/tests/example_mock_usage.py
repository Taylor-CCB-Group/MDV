#!/usr/bin/env python3
"""
Example usage of the mock AnnData module.

This script demonstrates how to use the MockAnnDataFactory and related utilities
for creating test data and stress testing the MDV conversion pipeline.

To create an MDV project with a unique column for manual integration testing
(edit from the frontend, then save state):

  python -m mdvtools.tests.example_mock_usage --create-unique-project ./test_unique_project

Then serve the project and open in the browser to edit the unique column and verify
round-trip / set_column_with_raw_data behaviour.
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


def example_create_project_with_unique_column(output_path=None):
    """Create an MDV project that has at least one unique column for manual integration testing.

    Use this to get a project you can serve and then edit unique columns from the frontend
    (e.g. change a cell_id value and save state) to verify unique column round-trip and
    set_column_with_raw_data behaviour.

    If output_path is None, uses a temporary directory (project is removed when the
    script exits). Pass a path (e.g. ./test_unique_project) to keep the project for
    serving and manual testing.

    To serve and test in the browser:
        cd /path/to/mdv && python -m mdvtools.serve --project /path/to/project
        # or use your usual serve command; then open the app and edit a unique column.
    """
    print("\n=== Create Project With Unique Column (for manual integration testing) ===")
    if output_path is None:
        project_dir = tempfile.mkdtemp(prefix="mdv_unique_")
        cleanup = True
    else:
        project_dir = os.path.expanduser(output_path)
        cleanup = False
        if os.path.exists(project_dir):
            shutil.rmtree(project_dir)
        os.makedirs(project_dir, exist_ok=True)

    try:
        factory = MockAnnDataFactory(random_seed=42)
        adata = factory.create_minimal(30, 15)
        with suppress_anndata_warnings():
            mdv = convert_scanpy_to_mdv(project_dir, adata, delete_existing=True)
        n_rows = len(adata)
        cell_ids = [f"cell_{i}" for i in range(n_rows)]
        mdv.set_column(
            "cells",
            {"name": "cell_id", "field": "cell_id", "datatype": "unique", "editable": True},
            cell_ids,
        )
        print(f"   Project created at: {os.path.abspath(project_dir)}")
        print(f"   Datasource 'cells' has a unique column 'cell_id' ({n_rows} rows).")
        if cleanup:
            print("   (Temporary directory; remove when done or re-run with output_path= to keep.)")
        else:
            print("   To serve for manual testing:")
            print(f"      python -m mdvtools.serve --project {os.path.abspath(project_dir)}")
    finally:
        if cleanup and os.path.exists(project_dir):
            shutil.rmtree(project_dir, ignore_errors=True)


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
            
            print("   Conversion successful")
            print(f"   Cells columns: {len(cells_metadata['columns'])}")
            print(f"   Genes columns: {len(genes_metadata['columns'])}")


def main():
    """Run all examples, or create a project with a unique column for manual testing."""
    import sys
    if len(sys.argv) >= 2 and sys.argv[1] == "--create-unique-project":
        output_path = sys.argv[2] if len(sys.argv) > 2 else "./test_unique_project"
        print("Creating MDV project with unique column for manual integration testing.")
        example_create_project_with_unique_column(output_path)
        return
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