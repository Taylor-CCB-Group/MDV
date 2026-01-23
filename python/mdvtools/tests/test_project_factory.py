#!/usr/bin/env python3
"""
Minimal test project factory for generating test data programmatically.

This module provides simple utilities to create test .h5ad files and MDV project
zips from either mock data or scanpy datasets, avoiding the need to commit large
binary test data files to the repository.
"""

import os
import tempfile
import logging
import shutil
import scanpy as sc
from mdvtools.tests.mock_anndata import create_minimal_anndata
from mdvtools.conversions import convert_scanpy_to_mdv

logger = logging.getLogger(__name__)


def create_test_h5ad_file(
    output_path: str,
    source: str = 'mock',
    dataset: str | None = None,
    n_cells: int = 100,
    n_genes: int = 200
) -> str:
    """Create a test .h5ad file for upload tests.
    
    Args:
        output_path: Where to save the file
        source: 'mock' or 'scanpy'
        dataset: For scanpy: 'pbmc3k' or 'pbmc3k_processed'
        n_cells: For mock data, number of cells
        n_genes: For mock data, number of genes
        
    Returns:
        Path to the created .h5ad file
        
    Example:
        >>> path = create_test_h5ad_file('/tmp/test.h5ad', source='mock', n_cells=100)
        >>> # Use in tests...
    """
    logger.info(f"Creating test h5ad file at {output_path} (source={source})")
    
    # Create AnnData
    if source == 'mock':
        adata = create_minimal_anndata(n_cells, n_genes)
    elif source == 'scanpy':
        if dataset == 'pbmc3k':
            adata = sc.datasets.pbmc3k()
        elif dataset == 'pbmc3k_processed':
            adata = sc.datasets.pbmc3k_processed()
        else:
            raise ValueError(f"Unknown dataset: {dataset}. Available: pbmc3k, pbmc3k_processed")
    else:
        raise ValueError(f"Unknown source: {source}. Available: mock, scanpy")
    
    # Ensure directory exists
    output_dir = os.path.dirname(output_path)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    
    # Write file
    adata.write_h5ad(output_path)
    logger.info(f"Created h5ad file with {adata.n_obs} cells and {adata.n_vars} genes")
    
    return output_path


def create_test_project_zip(
    output_path: str,
    source: str = 'mock',
    dataset: str | None = None,
    n_cells: int = 100,
    n_genes: int = 200
) -> str:
    """Create a test MDV project zip file for import tests.
    
    Args:
        output_path: Where to save the zip file
        source: 'mock' or 'scanpy'
        dataset: For scanpy: 'pbmc3k' or 'pbmc3k_processed'
        n_cells: For mock data, number of cells
        n_genes: For mock data, number of genes
        
    Returns:
        Path to the created .mdv.zip file
        
    Example:
        >>> path = create_test_project_zip('/tmp/test.mdv.zip', source='mock')
        >>> # Use in import tests...
    """
    logger.info(f"Creating test project zip at {output_path} (source={source})")
    
    # Create in temp dir
    with tempfile.TemporaryDirectory() as temp_dir:
        project_dir = os.path.join(temp_dir, 'project')
        
        # Create AnnData and convert to MDV
        if source == 'mock':
            adata = create_minimal_anndata(n_cells, n_genes)
        elif source == 'scanpy':
            if dataset == 'pbmc3k':
                adata = sc.datasets.pbmc3k()
            elif dataset == 'pbmc3k_processed':
                adata = sc.datasets.pbmc3k_processed()
            else:
                raise ValueError(f"Unknown dataset: {dataset}. Available: pbmc3k, pbmc3k_processed")
        else:
            raise ValueError(f"Unknown source: {source}. Available: mock, scanpy")
        
        logger.info("Converting AnnData to MDV project...")
        convert_scanpy_to_mdv(project_dir, adata)
        
        # Ensure output directory exists
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        
        # Create zip
        base_name = output_path.replace('.zip', '').replace('.mdv.zip', '')
        zip_path = shutil.make_archive(base_name, 'zip', project_dir)
        
        # Rename to .mdv.zip if needed
        if output_path.endswith('.mdv.zip') and not zip_path.endswith('.mdv.zip'):
            final_path = zip_path.replace('.zip', '.mdv.zip')
            shutil.move(zip_path, final_path)
            zip_path = final_path
        
        logger.info(f"Created project zip at {zip_path}")
    
    return zip_path
