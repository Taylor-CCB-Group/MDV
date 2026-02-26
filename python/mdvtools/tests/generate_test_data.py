#!/usr/bin/env python3
"""
Generate test MDV projects from mock or scanpy data.

Usage (with active mdv environment):
    python -m mdvtools.tests.generate_test_data ~/mdv/test_mock --mock
    python -m mdvtools.tests.generate_test_data ~/mdv/test_pbmc3k --scanpy pbmc3k_processed
    python -m mdvtools.tests.generate_test_data ~/mdv/test_large --mock --n-cells 1000000
"""

import argparse
import os
import sys
import inspect

import scanpy as sc
from mdvtools.tests.mock_anndata import create_minimal_anndata
from mdvtools.conversions import convert_scanpy_to_mdv


def _add_unique_column_to_project(project, n_rows):
    """Add a unique-typed column for manual testing of unique column edit/round-trip.
    Prefers datasource 'cells' when present so row count matches AnnData n_obs.
    """
    names = project.get_datasource_names()
    if not names:
        return
    ds_name = "cells" if "cells" in names else names[0]
    cell_ids = [f"cell_{i}" for i in range(n_rows)]
    project.set_column(
        ds_name,
        {"name": "cell_id", "field": "cell_id", "datatype": "unique", "editable": True},
        cell_ids,
    )


def main():
    parser = argparse.ArgumentParser(
        description='Generate test MDV projects from mock or scanpy data.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Create project from mock data (small, fast)
  python -m mdvtools.tests.generate_test_data ~/mdv/test_mock --mock

  # Create project from mock data (large)
  python -m mdvtools.tests.generate_test_data ~/mdv/test_large --mock --n-cells 1000000 --n-genes 5000

  # Create project with a unique column for manual integration testing (edit from frontend)
  python -m mdvtools.tests.generate_test_data ~/mdv/test_unique --mock --with-unique-column

  # Create project from scanpy dataset
  python -m mdvtools.tests.generate_test_data ~/mdv/test_pbmc3k --scanpy pbmc3k_processed

See https://scanpy.readthedocs.io/en/1.11.x/api/datasets.html for available scanpy datasets.
        """
    )
    
    parser.add_argument(
        'output_path',
        type=str,
        help='Output path for the MDV project (e.g., ~/mdv/test_project)'
    )
    
    # Data source (mutually exclusive)
    source_group = parser.add_mutually_exclusive_group(required=True)
    source_group.add_argument(
        '--mock',
        action='store_true',
        help='Generate from mock AnnData'
    )
    source_group.add_argument(
        '--scanpy',
        type=str,
        metavar='DATASET',
        help='Load from scanpy.datasets (any available dataset)'
    )
    
    # Mock data parameters
    parser.add_argument(
        '--n-cells',
        type=int,
        default=100,
        help='Number of cells for mock data (default: 100)'
    )
    parser.add_argument(
        '--n-genes',
        type=int,
        default=200,
        help='Number of genes for mock data (default: 200)'
    )
    
    # Other options
    parser.add_argument(
        '--force',
        action='store_true',
        help='Overwrite existing project'
    )
    parser.add_argument(
        '--with-unique-column',
        action='store_true',
        help='Add a unique-typed column (e.g. cell_id) to the first datasource for manual testing of unique column edit/round-trip'
    )
    
    args = parser.parse_args()
    
    # Expand user path
    output_path = os.path.expanduser(args.output_path)
    
    # Check if exists
    if os.path.exists(output_path) and not args.force:
        print(f"Error: {output_path} already exists. Use --force to overwrite.")
        sys.exit(1)
    
    # Create AnnData from source
    if args.mock:
        print(f"Creating mock AnnData with {args.n_cells:,} cells and {args.n_genes:,} genes...")
        adata = create_minimal_anndata(n_cells=args.n_cells, n_genes=args.n_genes)
    else:
        dataset = args.scanpy
        print(f"Loading scanpy dataset: {dataset}...")
        
        # Dynamically get all available datasets from sc.datasets
        # notwithstanding ebi_expression_atlas or anything else which might need any arguments
        dataset_loaders = {
            name: loader
            for name in dir(sc.datasets)
            if not name.startswith('_')
            and callable(loader := getattr(sc.datasets, name, None))
            and not any(
                p.default is p.empty and p.kind in (p.POSITIONAL_ONLY, p.POSITIONAL_OR_KEYWORD)
                for p in inspect.signature(loader).parameters.values()
            )
        }
        if dataset not in dataset_loaders:
            available = ', '.join(sorted(dataset_loaders.keys()))
            print(f"Error: Unknown dataset '{dataset}'. Available: {available}")
            sys.exit(1)
        
        adata = dataset_loaders[dataset]()
        assert(isinstance(adata, sc.AnnData))
    
    print(f"AnnData: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    
    # Convert to MDV project
    print(f"Creating MDV project at {output_path}...")
    convert_scanpy_to_mdv(output_path, adata, delete_existing=args.force)
    
    # Add provenance metadata to state.json
    print("Adding provenance metadata...")
    from mdvtools.mdvproject import MDVProject
    from datetime import datetime
    
    project = MDVProject(output_path)
    state = project.state
    
    # Build provenance metadata
    provenance = {
        "created_by": "generate_test_data.py",
        "created_at": datetime.now().isoformat(),
        "source": "mock" if args.mock else "scanpy",
        "parameters": {}
    }
    
    # Add source-specific parameters
    if args.mock:
        provenance["parameters"]["n_cells"] = args.n_cells
        provenance["parameters"]["n_genes"] = args.n_genes
    else:
        provenance["parameters"]["dataset"] = args.scanpy
        # Also include actual dimensions from the loaded data
        provenance["parameters"]["n_cells"] = adata.n_obs
        provenance["parameters"]["n_genes"] = adata.n_vars
    
    # Add provenance to state
    state["provenance"] = provenance
    project.state = state
    
    if args.with_unique_column:
        n_rows = args.n_cells if args.mock else adata.n_obs
        _add_unique_column_to_project(project, n_rows)
        print("  Added unique column to first datasource for manual integration testing.")
    
    print(f"✓ Created MDV project at {output_path}")
    print(f"  Provenance: {provenance['source']} data, {provenance['parameters']}")


if __name__ == '__main__':
    main()
