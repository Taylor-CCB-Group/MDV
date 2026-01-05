"""Manages h5ad file for MCP-Bio integration.

This module handles the lifecycle of h5ad files used by MCP-Bio tools,
including finding existing files, creating them from MDV datasources,
and syncing results back to MDV.
"""
import os
import logging
from typing import List, Optional

import numpy as np

logger = logging.getLogger(__name__)

# Lazy imports to avoid hard dependencies
def _import_scanpy():
    try:
        import scanpy as sc
        return sc
    except ImportError:
        raise ImportError(
            "scanpy is required for MCP-Bio integration. "
            "Install it with: pip install scanpy"
        )

def _import_anndata():
    try:
        from anndata import AnnData
        return AnnData
    except ImportError:
        raise ImportError(
            "anndata is required for MCP-Bio integration. "
            "Install it with: pip install anndata"
        )


class H5ADManager:
    """Manages h5ad file lifecycle for MCP-Bio tools.
    
    This class handles:
    - Finding existing h5ad files in the MDV project
    - Creating h5ad files from MDV datasources
    - Syncing analysis results back to MDV
    """
    
    def __init__(self, project):
        """Initialize the H5AD manager.
        
        Args:
            project: MDVProject instance
        """
        from mdvtools.mdvproject import MDVProject
        if not isinstance(project, MDVProject):
            raise TypeError(f"Expected MDVProject, got {type(project)}")
        
        self.project = project
        self._h5ad_path: Optional[str] = None
        self._initialized = False
    
    def _find_or_create_h5ad(self) -> str:
        """Find existing h5ad or create from MDV datasources.
        
        Returns:
            Path to h5ad file
        """
        # First, look for existing h5ad files
        for file in os.listdir(self.project.dir):
            if file.endswith(".h5ad"):
                path = os.path.join(self.project.dir, file)
                logger.info(f"Found existing h5ad file: {path}")
                return path
        
        # No existing file, create from MDV
        logger.info("No h5ad file found, creating from MDV datasources")
        return self._create_from_mdv()
    
    def _create_from_mdv(self) -> str:
        """Create minimal h5ad from MDV datasources.
        
        Creates an AnnData object with:
        - obs from 'cells' datasource
        - var from 'genes' datasource (if available)
        - X matrix from expression data (if available)
        
        Returns:
            Path to newly created h5ad file
        """
        AnnData = _import_anndata()
        
        path = os.path.join(self.project.dir, "anndata.h5ad")
        
        # Get cells datasource
        cells_df = None
        genes_df = None
        
        for ds in self.project.datasources:
            ds_name = ds['name']
            if ds_name.lower() in ('cells', 'cell', 'obs'):
                cells_df = self.project.get_datasource_as_dataframe(ds_name)
                logger.info(f"Using '{ds_name}' as cells/obs")
            elif ds_name.lower() in ('genes', 'gene', 'var', 'features'):
                genes_df = self.project.get_datasource_as_dataframe(ds_name)
                logger.info(f"Using '{ds_name}' as genes/var")
        
        # Fallback to first datasource if no cells found
        if cells_df is None and len(self.project.datasources) > 0:
            ds_name = self.project.datasources[0]['name']
            cells_df = self.project.get_datasource_as_dataframe(ds_name)
            logger.warning(f"No 'cells' datasource found, using '{ds_name}'")
        
        if cells_df is None:
            raise ValueError("No datasources available to create h5ad file")
        
        # Create AnnData object
        if genes_df is not None:
            # Create with proper dimensions for single-cell data
            n_obs = len(cells_df)
            n_vars = len(genes_df)
            
            # Check if there's expression data we can include
            # For now, create empty matrix - actual expression data
            # should be loaded separately if needed
            X = np.zeros((n_obs, n_vars), dtype=np.float32)
            
            adata = AnnData(
                X=X,
                obs=cells_df.reset_index(drop=True),
                var=genes_df.reset_index(drop=True)
            )
        else:
            # Just obs data
            adata = AnnData(obs=cells_df.reset_index(drop=True))
        
        # Write to disk
        adata.write_h5ad(path)
        logger.info(f"Created h5ad file: {path} (obs: {adata.n_obs}, var: {adata.n_vars})")
        
        return path
    
    def get_path(self) -> str:
        """Get path to h5ad file, creating if necessary.
        
        Returns:
            Absolute path to h5ad file
        """
        if not self._initialized:
            self._h5ad_path = self._find_or_create_h5ad()
            self._initialized = True
        # #region agent log
        log_debug("H6", "H5AD path retrieved", {"path": self._h5ad_path})
        # #endregion
        return self._h5ad_path
    
    def sync_to_mdv(self, result_path: str, columns: List[str]):
        """Sync columns from result h5ad back to MDV.
        
        Args:
            result_path: Path to h5ad file with analysis results
            columns: List of column names to sync
        """
        sc = _import_scanpy()
        
        logger.info(f"Syncing columns from {result_path}: {columns}")
        
        # Read the result file
        adata = sc.read_h5ad(result_path)
        
        synced = []
        failed = []
        
        for col in columns:
            try:
                if col in adata.obs.columns:
                    # Regular obs column (e.g., 'leiden', 'cell_type')
                    data = adata.obs[col].values
                    
                    # Handle categorical data
                    if hasattr(data, 'codes'):
                        # Convert categorical to string for MDV
                        data = data.astype(str)
                    
                    self.project.set_column("cells", col, data)
                    synced.append(col)
                    logger.info(f"Synced obs column: {col}")
                    
                elif col.startswith(("umap_", "pca_", "tsne_")):
                    # Embedding coordinates from obsm
                    # Format: "umap_0", "umap_1", "pca_0", etc.
                    parts = col.rsplit("_", 1)
                    if len(parts) == 2:
                        embed = parts[0]  # e.g., "umap"
                        idx = int(parts[1])  # e.g., 0
                        key = f"X_{embed}"
                        
                        if key in adata.obsm:
                            data = adata.obsm[key][:, idx]
                            self.project.set_column("cells", col, data)
                            synced.append(col)
                            logger.info(f"Synced embedding column: {col} from {key}")
                        else:
                            logger.warning(f"Embedding {key} not found in obsm")
                            failed.append(col)
                    else:
                        logger.warning(f"Could not parse embedding column: {col}")
                        failed.append(col)
                else:
                    logger.warning(f"Column {col} not found in obs or recognized as embedding")
                    failed.append(col)
                    
            except Exception as e:
                logger.error(f"Failed to sync column {col}: {e}")
                failed.append(col)
        
        if failed:
            logger.warning(f"Failed to sync columns: {failed}")
        
        return {"synced": synced, "failed": failed}
    
    def refresh(self):
        """Force re-detection of h5ad file on next get_path() call."""
        self._initialized = False
        self._h5ad_path = None




# #region agent log
import json
import time
def log_debug(hypothesis_id, message, data=None):
    log_entry = {
        "id": f"log_{int(time.time())}_{hypothesis_id}",
        "timestamp": int(time.time() * 1000),
        "location": "h5ad_manager.py",
        "message": message,
        "data": data or {},
        "sessionId": "debug-session",
        "runId": "run1",
        "hypothesisId": hypothesis_id
    }
    try:
        with open("/app/.cursor/debug.log", "a") as f:
            f.write(json.dumps(log_entry) + "\n")
    except Exception:
        pass
# #endregion
