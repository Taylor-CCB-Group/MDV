#!/usr/bin/env python3
"""
Integrate external metadata (CSV/TSV) with MDV projects.

Usage:
    python -m mdvtools.umts.integrate_metadata \
        --project /app/mdv/my_project \
        --metadata cell_phenotypes.csv \
        --key cell_id \
        --datasource cells
"""

import argparse
import pandas as pd
from pathlib import Path
import sys

def integrate_metadata(
    project_path: str,
    metadata_file: str,
    join_key: str = "cell_id",
    metadata_key: str = None,
    datasource: str = "cells",
    overwrite: bool = False
):
    """
    Integrate external metadata into an MDV project.
    
    Args:
        project_path: Path to MDV project directory
        metadata_file: Path to CSV/TSV file with metadata
        join_key: Column name in MDV datasource to join on
        metadata_key: Column name in metadata file to join on (defaults to join_key)
        datasource: MDV datasource name (default: "cells")
        overwrite: Whether to overwrite existing columns
    """
    from mdvtools import MDVProject
    
    # Detect separator
    sep = '\t' if metadata_file.endswith('.tsv') else ','
    
    print(f"Loading metadata from {metadata_file}...")
    metadata = pd.read_csv(metadata_file, sep=sep)
    
    if metadata_key is None:
        metadata_key = join_key
    
    print(f"Loading MDV project from {project_path}...")
    mdv = MDVProject(project_path)
    
    # Get join column from project
    print(f"Getting {join_key} from {datasource}...")
    project_keys = mdv.get_column(datasource, join_key)
    
    # Check for join key in metadata
    if metadata_key not in metadata.columns:
        raise ValueError(f"Column '{metadata_key}' not found in metadata file")
    
    # Create lookup dictionary for each metadata column
    metadata_dict = metadata.set_index(metadata_key)
    
    # Get columns to add (exclude the key column)
    columns_to_add = [col for col in metadata.columns if col != metadata_key]
    
    print(f"\nAdding {len(columns_to_add)} metadata columns:")
    for col in columns_to_add:
        print(f"  - {col}")
    
    # Check for existing columns
    existing_cols = [col for col in columns_to_add 
                     if any(c['field'] == col for c in mdv.get_datasource_metadata(datasource)['columns'])]
    
    if existing_cols and not overwrite:
        print(f"\n⚠ Warning: These columns already exist: {existing_cols}")
        print("Use --overwrite to replace them, or they will be skipped.")
        columns_to_add = [col for col in columns_to_add if col not in existing_cols]
    
    # Add each column
    added_count = 0
    for col in columns_to_add:
        # Match values to project keys
        values = [metadata_dict.loc[key, col] if key in metadata_dict.index 
                  else None for key in project_keys]
        
        # Add column to project
        mdv.set_column(datasource, col, values)
        added_count += 1
        print(f"  ✓ Added {col}")
    
    print(f"\n✓ Successfully added {added_count} columns to {datasource}")
    
    # Show summary
    matched = sum(1 for key in project_keys if key in metadata_dict.index)
    print(f"  Matched: {matched}/{len(project_keys)} rows")
    
    return mdv


def main():
    parser = argparse.ArgumentParser(
        description='Integrate external metadata into MDV projects'
    )
    
    parser.add_argument('--project', required=True,
                       help='Path to MDV project directory')
    parser.add_argument('--metadata', required=True,
                       help='Path to metadata CSV/TSV file')
    parser.add_argument('--key', default='cell_id',
                       help='Column name in MDV datasource to join on (default: cell_id)')
    parser.add_argument('--metadata-key', default=None,
                       help='Column name in metadata file to join on (default: same as --key)')
    parser.add_argument('--datasource', default='cells',
                       help='MDV datasource name (default: cells)')
    parser.add_argument('--overwrite', action='store_true',
                       help='Overwrite existing columns')
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='Verbose output')
    
    args = parser.parse_args()
    
    # Validate paths
    if not Path(args.project).exists():
        print(f"❌ Error: Project not found: {args.project}")
        sys.exit(1)
    
    if not Path(args.metadata).exists():
        print(f"❌ Error: Metadata file not found: {args.metadata}")
        sys.exit(1)
    
    try:
        integrate_metadata(
            project_path=args.project,
            metadata_file=args.metadata,
            join_key=args.key,
            metadata_key=args.metadata_key,
            datasource=args.datasource,
            overwrite=args.overwrite
        )
        
        print("\n✓ Metadata integration complete!")
        return 0
        
    except Exception as e:
        print(f"\n❌ Error: {str(e)}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())

