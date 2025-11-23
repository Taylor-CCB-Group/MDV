#!/usr/bin/env python3
"""
UMTS Command Line Interface for MDVTools.

Provides command-line tools for converting various data formats to MDV projects
with automatic lineage tracking.
"""

import argparse
import sys
import os
from pathlib import Path


def cmd_convert_seurat(args):
    """Convert Seurat RDS file to MDV project."""
    from .seurat_reader import convert_seurat_to_mdv, check_seurat_dependencies
    
    # Check dependencies first
    print("Checking dependencies...")
    deps = check_seurat_dependencies()
    
    if not deps['ready']:
        print("\n❌ Missing required dependencies:")
        if not deps['rpy2']:
            print("  - rpy2 (install: pip install rpy2)")
        if not deps['anndata2ri']:
            print("  - anndata2ri (install: pip install anndata2ri)")
        if not deps['r_available']:
            print("  - R (install R from https://www.r-project.org/)")
        if not deps['seurat_available']:
            print("  - Seurat R package (install in R: install.packages('Seurat'))")
        if not deps['seuratdisk_available']:
            print("  - SeuratDisk R package (install in R: remotes::install_github('mojaveazure/seurat-disk'))")
        print("\nPlease install missing dependencies and try again.")
        sys.exit(1)
    
    print("✓ All dependencies available\n")
    
    # Convert
    try:
        mdv = convert_seurat_to_mdv(
            seurat_file=args.input,
            output_folder=args.output,
            max_dims=args.max_dims,
            delete_existing=not args.no_delete,
            label=args.label,
            chunk_data=args.chunk_data,
            track_lineage=not args.no_lineage
        )
        
        print(f"\n✓ Successfully created MDV project at: {args.output}")
        
        if not args.no_lineage:
            lineage = mdv.get_lineage()
            if lineage:
                print(f"✓ Lineage information saved")
                print(f"  - Source: {lineage['source_files'][0]['path']}")
                print(f"  - SHA256: {lineage['source_files'][0]['sha256'][:16]}...")
                print(f"  - Timestamp: {lineage['created_timestamp']}")
        
        return 0
        
    except Exception as e:
        print(f"\n❌ Error: {str(e)}", file=sys.stderr)
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1


def cmd_convert_h5ad(args):
    """Convert h5ad file to MDV project."""
    import scanpy as sc
    from ..conversions import convert_scanpy_to_mdv
    
    try:
        print(f"Loading AnnData from {args.input}...")
        adata = sc.read_h5ad(args.input)
        print(f"Loaded: {adata.n_obs} cells, {adata.n_vars} genes\n")
        
        print(f"Converting to MDV project at {args.output}...")
        mdv = convert_scanpy_to_mdv(
            folder=args.output,
            scanpy_object=adata,
            max_dims=args.max_dims,
            delete_existing=not args.no_delete,
            label=args.label,
            chunk_data=args.chunk_data,
            track_lineage=not args.no_lineage,
            source_file=args.input
        )
        
        print(f"\n✓ Successfully created MDV project at: {args.output}")
        
        if not args.no_lineage:
            lineage = mdv.get_lineage()
            if lineage:
                print(f"✓ Lineage information saved")
                print(f"  - Source: {lineage['source_files'][0]['path']}")
                print(f"  - SHA256: {lineage['source_files'][0]['sha256'][:16]}...")
                print(f"  - Timestamp: {lineage['created_timestamp']}")
        
        return 0
        
    except Exception as e:
        print(f"\n❌ Error: {str(e)}", file=sys.stderr)
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1


def cmd_show_lineage(args):
    """Show lineage information for an MDV project."""
    from ..mdvproject import MDVProject
    
    try:
        mdv = MDVProject(args.project)
        lineage = mdv.get_lineage()
        
        if lineage is None:
            print(f"No lineage information found in {args.project}")
            return 1
        
        print("=" * 60)
        print(f"LINEAGE INFORMATION: {args.project}")
        print("=" * 60)
        
        print(f"\nUMTS Version: {lineage['umts_version']}")
        print(f"Created: {lineage['created_timestamp']}")
        
        if lineage.get('source_files'):
            print(f"\nSource Files ({len(lineage['source_files'])}):")
            for i, src in enumerate(lineage['source_files'], 1):
                print(f"  {i}. {src['path']}")
                print(f"     SHA256: {src['sha256']}")
                print(f"     Size: {src['size_bytes']:,} bytes")
                print(f"     Modified: {src['modified_timestamp']}")
        
        print(f"\nConversion:")
        print(f"  Function: {lineage['conversion']['function']}")
        print(f"  Parameters:")
        for key, value in lineage['conversion']['parameters'].items():
            print(f"    {key}: {value}")
        
        print(f"\nEnvironment:")
        print(f"  Python: {lineage['environment']['python_version']}")
        print(f"  Platform: {lineage['environment']['platform']}")
        
        if args.packages:
            print(f"\n  Packages:")
            for pkg, ver in lineage['environment']['packages'].items():
                print(f"    {pkg}: {ver}")
        
        if lineage.get('notes'):
            print(f"\nNotes ({len(lineage['notes'])}):")
            for note in lineage['notes']:
                print(f"  [{note['timestamp']}] {note['message']}")
        
        print("\n" + "=" * 60)
        return 0
        
    except Exception as e:
        print(f"❌ Error: {str(e)}", file=sys.stderr)
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1


def cmd_check_deps(args):
    """Check Seurat conversion dependencies."""
    from .seurat_reader import check_seurat_dependencies
    
    print("Checking UMTS dependencies...\n")
    
    deps = check_seurat_dependencies()
    
    def status_icon(available):
        return "✓" if available else "✗"
    
    print("Python Packages:")
    print(f"  {status_icon(deps['rpy2'])} rpy2")
    print(f"  {status_icon(deps['anndata2ri'])} anndata2ri")
    
    print("\nR Environment:")
    print(f"  {status_icon(deps['r_available'])} R available")
    print(f"  {status_icon(deps['seurat_available'])} Seurat package")
    print(f"  {status_icon(deps['seuratdisk_available'])} SeuratDisk package")
    
    print(f"\nOverall Status: {'✓ Ready' if deps['ready'] else '✗ Not Ready'}")
    
    if not deps['ready']:
        print("\nTo enable Seurat conversion, install missing dependencies:")
        if not deps['rpy2'] or not deps['anndata2ri']:
            print("  pip install rpy2 anndata2ri")
        if not deps['r_available']:
            print("  Install R from: https://www.r-project.org/")
        if not deps['seurat_available']:
            print("  In R: install.packages('Seurat')")
        if not deps['seuratdisk_available']:
            print("  In R: remotes::install_github('mojaveazure/seurat-disk')")
        return 1
    
    return 0


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        prog='umts',
        description='UMTS - Universal Modality Translator System for MDVTools',
        epilog='For more information, see: docs/UMTS_LINEAGE.md'
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # Convert Seurat command
    parser_seurat = subparsers.add_parser(
        'convert-seurat',
        help='Convert Seurat RDS file to MDV project',
        description='Convert a Seurat object from an RDS file to MDV format with lineage tracking'
    )
    parser_seurat.add_argument('input', help='Input Seurat RDS file')
    parser_seurat.add_argument('output', help='Output MDV project directory')
    parser_seurat.add_argument('--max-dims', type=int, default=3,
                              help='Maximum dimensions to include (default: 3)')
    parser_seurat.add_argument('--label', default='',
                              help='Prefix for datasource names (default: none)')
    parser_seurat.add_argument('--chunk-data', action='store_true',
                              help='Process data in chunks (slower, lower memory)')
    parser_seurat.add_argument('--no-delete', action='store_true',
                              help='Do not delete existing project (merge instead)')
    parser_seurat.add_argument('--no-lineage', action='store_true',
                              help='Disable lineage tracking')
    parser_seurat.add_argument('-v', '--verbose', action='store_true',
                              help='Verbose output with stack traces')
    parser_seurat.set_defaults(func=cmd_convert_seurat)
    
    # Convert h5ad command
    parser_h5ad = subparsers.add_parser(
        'convert-h5ad',
        help='Convert h5ad file to MDV project',
        description='Convert an AnnData h5ad file to MDV format with lineage tracking'
    )
    parser_h5ad.add_argument('input', help='Input h5ad file')
    parser_h5ad.add_argument('output', help='Output MDV project directory')
    parser_h5ad.add_argument('--max-dims', type=int, default=3,
                            help='Maximum dimensions to include (default: 3)')
    parser_h5ad.add_argument('--label', default='',
                            help='Prefix for datasource names (default: none)')
    parser_h5ad.add_argument('--chunk-data', action='store_true',
                            help='Process data in chunks (slower, lower memory)')
    parser_h5ad.add_argument('--no-delete', action='store_true',
                            help='Do not delete existing project (merge instead)')
    parser_h5ad.add_argument('--no-lineage', action='store_true',
                            help='Disable lineage tracking')
    parser_h5ad.add_argument('-v', '--verbose', action='store_true',
                            help='Verbose output with stack traces')
    parser_h5ad.set_defaults(func=cmd_convert_h5ad)
    
    # Show lineage command
    parser_lineage = subparsers.add_parser(
        'show-lineage',
        help='Show lineage information for an MDV project',
        description='Display provenance and lineage information for an MDV project'
    )
    parser_lineage.add_argument('project', help='MDV project directory')
    parser_lineage.add_argument('--packages', action='store_true',
                               help='Show all package versions')
    parser_lineage.add_argument('-v', '--verbose', action='store_true',
                               help='Verbose output with stack traces')
    parser_lineage.set_defaults(func=cmd_show_lineage)
    
    # Check dependencies command
    parser_deps = subparsers.add_parser(
        'check-deps',
        help='Check Seurat conversion dependencies',
        description='Verify that all required dependencies for Seurat conversion are installed'
    )
    parser_deps.set_defaults(func=cmd_check_deps)
    
    # Integrate metadata command
    parser_metadata = subparsers.add_parser(
        'add-metadata',
        help='Add external metadata to existing MDV project',
        description='Integrate CSV/TSV metadata into an existing MDV project'
    )
    parser_metadata.add_argument('project', help='MDV project directory')
    parser_metadata.add_argument('metadata', help='CSV/TSV metadata file')
    parser_metadata.add_argument('--key', default='cell_id',
                                help='Join key in MDV datasource (default: cell_id)')
    parser_metadata.add_argument('--metadata-key', 
                                help='Join key in metadata file (default: same as --key)')
    parser_metadata.add_argument('--datasource', default='cells',
                                help='MDV datasource name (default: cells)')
    parser_metadata.add_argument('--overwrite', action='store_true',
                                help='Overwrite existing columns')
    parser_metadata.add_argument('-v', '--verbose', action='store_true',
                                help='Verbose output')
    
    def cmd_add_metadata(args):
        """Add metadata to existing project."""
        from .integrate_metadata import integrate_metadata
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
            print(f"\n❌ Error: {str(e)}", file=sys.stderr)
            if args.verbose:
                import traceback
                traceback.print_exc()
            return 1
    
    parser_metadata.set_defaults(func=cmd_add_metadata)
    
    # Parse arguments
    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        return 1
    
    # Execute command
    return args.func(args)


if __name__ == '__main__':
    sys.exit(main())

