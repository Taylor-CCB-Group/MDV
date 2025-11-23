"""
Seurat RDS file reader for UMTS.

Provides functionality to read Seurat objects from R and convert them to AnnData
for use with MDVTools.
"""

import os
from typing import Optional


def convert_seurat_to_mdv(
    seurat_file: str,
    output_folder: str,
    max_dims: int = 3,
    delete_existing: bool = True,
    label: str = "",
    chunk_data: bool = False,
    add_layer_data: bool = True,
    gene_identifier_column: Optional[str] = None,
    track_lineage: bool = True
):
    """
    Convert a Seurat RDS file directly to MDV format.
    
    This function reads a Seurat object from an RDS file, converts it to AnnData,
    and then creates an MDV project with lineage tracking.
    
    Args:
        seurat_file: Path to the Seurat RDS file
        output_folder: Path to output MDV project folder
        max_dims: Maximum number of dimensions to include from dimensionality reductions
        delete_existing: Whether to delete existing project data
        label: Prefix to add to datasource names
        chunk_data: Whether to process data in chunks (slower but lower memory)
        add_layer_data: Whether to include layer data
        gene_identifier_column: Column to use for gene identifiers
        track_lineage: Whether to track lineage (default: True)
        
    Returns:
        MDVProject: The created MDV project
        
    Raises:
        ImportError: If required packages (rpy2, anndata2ri) are not installed
        FileNotFoundError: If the Seurat file doesn't exist
        RuntimeError: If conversion fails
    """
    # Check file exists
    if not os.path.exists(seurat_file):
        raise FileNotFoundError(f"Seurat file not found: {seurat_file}")
    
    # Try to import required packages
    try:
        import rpy2.robjects as ro
        from rpy2.robjects import pandas2ri, numpy2ri
        from rpy2.robjects.conversion import localconverter
        import anndata2ri
    except ImportError as e:
        raise ImportError(
            "Seurat RDS reading requires rpy2 and anndata2ri packages. "
            "Install with: pip install rpy2 anndata2ri\n"
            "You also need R installed with Seurat and SeuratDisk packages."
        ) from e
    
    try:
        # Activate converters
        pandas2ri.activate()
        numpy2ri.activate()
        anndata2ri.activate()
        
        print(f"Loading Seurat object from {seurat_file}...")
        
        # Load R libraries
        ro.r('library(Seurat)')
        ro.r('library(SeuratDisk)')
        
        # Read the RDS file
        ro.r(f'seurat_obj <- readRDS("{seurat_file}")')
        
        # Create temporary h5Seurat file
        import tempfile
        with tempfile.TemporaryDirectory() as tmpdir:
            h5seurat_file = os.path.join(tmpdir, "temp.h5Seurat")
            h5ad_file = os.path.join(tmpdir, "temp.h5ad")
            
            print("Converting Seurat to h5ad format...")
            
            # Save as h5Seurat
            ro.r(f'SaveH5Seurat(seurat_obj, filename = "{h5seurat_file}", overwrite = TRUE)')
            
            # Convert to h5ad
            ro.r(f'Convert("{h5seurat_file}", dest = "h5ad", overwrite = TRUE)')
            
            # Load as AnnData
            print("Loading as AnnData...")
            import scanpy as sc
            adata = sc.read_h5ad(h5ad_file)
            
            print(f"Loaded: {adata.n_obs} cells, {adata.n_vars} genes")
            
            # Convert to MDV with lineage tracking
            from ..conversions import convert_scanpy_to_mdv
            
            print(f"Converting to MDV project at {output_folder}...")
            mdv = convert_scanpy_to_mdv(
                folder=output_folder,
                scanpy_object=adata,
                max_dims=max_dims,
                delete_existing=delete_existing,
                label=label,
                chunk_data=chunk_data,
                add_layer_data=add_layer_data,
                gene_identifier_column=gene_identifier_column,
                track_lineage=track_lineage,
                source_file=seurat_file  # Track original Seurat file
            )
            
            print(f"✓ Successfully created MDV project")
            
            return mdv
            
    except Exception as e:
        raise RuntimeError(f"Error converting Seurat file: {str(e)}") from e
    finally:
        # Deactivate converters
        try:
            pandas2ri.deactivate()
            numpy2ri.deactivate()
            anndata2ri.deactivate()
        except:
            pass


def check_seurat_dependencies() -> dict:
    """
    Check if required dependencies for Seurat reading are available.
    
    Returns:
        Dictionary with dependency status:
        {
            'rpy2': bool,
            'anndata2ri': bool,
            'r_available': bool,
            'seurat_available': bool,
            'seuratdisk_available': bool,
            'ready': bool  # True if all dependencies available
        }
    """
    status = {
        'rpy2': False,
        'anndata2ri': False,
        'r_available': False,
        'seurat_available': False,
        'seuratdisk_available': False,
        'ready': False
    }
    
    # Check Python packages
    try:
        import rpy2
        status['rpy2'] = True
    except ImportError:
        pass
    
    try:
        import anndata2ri
        status['anndata2ri'] = True
    except ImportError:
        pass
    
    # Check R availability
    if status['rpy2']:
        try:
            import rpy2.robjects as ro
            ro.r('R.version')
            status['r_available'] = True
            
            # Check for Seurat
            try:
                ro.r('library(Seurat)')
                status['seurat_available'] = True
            except:
                pass
            
            # Check for SeuratDisk
            try:
                ro.r('library(SeuratDisk)')
                status['seuratdisk_available'] = True
            except:
                pass
                
        except:
            pass
    
    status['ready'] = all([
        status['rpy2'],
        status['anndata2ri'],
        status['r_available'],
        status['seurat_available'],
        status['seuratdisk_available']
    ])
    
    return status

