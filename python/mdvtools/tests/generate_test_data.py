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

import scanpy as sc
from mdvtools.tests.mock_anndata import create_minimal_anndata
from mdvtools.conversions import convert_scanpy_to_mdv


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
        help='Load from scanpy.datasets (pbmc3k, pbmc3k_processed, pbmc68k_reduced)'
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
        
        dataset_loaders = {
            'pbmc3k': sc.datasets.pbmc3k,
            'pbmc3k_processed': sc.datasets.pbmc3k_processed,
            'pbmc68k_reduced': sc.datasets.pbmc68k_reduced,
        }
        
        if dataset not in dataset_loaders:
            available = ', '.join(dataset_loaders.keys())
            print(f"Error: Unknown dataset '{dataset}'. Available: {available}")
            sys.exit(1)
        
        adata = dataset_loaders[dataset]()
    
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
    
    print(f"✓ Created MDV project at {output_path}")
    print(f"  Provenance: {provenance['source']} data, {provenance['parameters']}")


if __name__ == '__main__':
    main()
