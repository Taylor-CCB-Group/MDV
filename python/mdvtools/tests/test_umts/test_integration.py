"""
Integration tests for UMTS lineage tracking with MDV conversions.
"""

import pytest
import os
import json
import tempfile
import numpy as np
import pandas as pd
from pathlib import Path

try:
    from scanpy import AnnData
    SCANPY_AVAILABLE = True
except ImportError:
    SCANPY_AVAILABLE = False

from mdvtools.conversions import convert_scanpy_to_mdv
from mdvtools.mdvproject import MDVProject


@pytest.mark.skipif(not SCANPY_AVAILABLE, reason="scanpy not installed")
class TestConversionWithLineage:
    """Test convert_scanpy_to_mdv with lineage tracking."""
    
    def create_mock_anndata(self, n_obs=100, n_vars=50):
        """Create a minimal mock AnnData object for testing."""
        # Create mock data
        X = np.random.rand(n_obs, n_vars)
        
        obs = pd.DataFrame({
            'cell_type': [f'type_{i % 5}' for i in range(n_obs)],
            'batch': [f'batch_{i % 2}' for i in range(n_obs)]
        }, index=[f'cell_{i}' for i in range(n_obs)])
        
        var = pd.DataFrame({
            'gene_name': [f'gene_{i}' for i in range(n_vars)],
            'highly_variable': np.random.choice([True, False], n_vars)
        }, index=[f'GENE{i}' for i in range(n_vars)])
        
        # Create AnnData
        adata = AnnData(X=X, obs=obs, var=var)
        
        # Add a dimensionality reduction (UMAP)
        adata.obsm['X_umap'] = np.random.rand(n_obs, 2)
        
        return adata
    
    def test_conversion_with_lineage_enabled(self, tmp_path):
        """Test that lineage tracking works when enabled."""
        # Create mock data
        adata = self.create_mock_anndata()
        
        # Save to file for source tracking
        source_file = tmp_path / "test.h5ad"
        adata.write_h5ad(str(source_file))
        
        # Convert with lineage tracking
        project_dir = tmp_path / "test_project"
        mdv = convert_scanpy_to_mdv(
            folder=str(project_dir),
            scanpy_object=adata,
            max_dims=2,
            delete_existing=True,
            track_lineage=True,
            source_file=str(source_file)
        )
        
        # Check that lineage file was created
        lineage_file = project_dir / "lineage.json"
        assert lineage_file.exists(), "lineage.json should be created"
        
        # Load and verify lineage data
        with open(lineage_file) as f:
            lineage = json.load(f)
        
        assert 'umts_version' in lineage
        assert 'created_timestamp' in lineage
        assert 'source_files' in lineage
        assert 'conversion' in lineage
        assert 'environment' in lineage
        
        # Check source file info
        assert len(lineage['source_files']) == 1
        source = lineage['source_files'][0]
        assert 'sha256' in source
        assert 'size_bytes' in source
        assert source['path'] == str(source_file.absolute())
        
        # Check conversion parameters
        params = lineage['conversion']['parameters']
        assert params['max_dims'] == 2
        assert params['delete_existing'] == True
        assert params['n_obs'] == 100
        assert params['n_vars'] == 50
        
        # Check environment
        assert 'python_version' in lineage['environment']
        assert 'packages' in lineage['environment']
        
    def test_conversion_with_lineage_disabled(self, tmp_path):
        """Test that lineage tracking can be disabled."""
        adata = self.create_mock_anndata(n_obs=50, n_vars=30)
        
        project_dir = tmp_path / "test_project_no_lineage"
        mdv = convert_scanpy_to_mdv(
            folder=str(project_dir),
            scanpy_object=adata,
            track_lineage=False
        )
        
        # Check that lineage file was NOT created
        lineage_file = project_dir / "lineage.json"
        assert not lineage_file.exists(), "lineage.json should not be created when disabled"
    
    def test_conversion_without_source_file(self, tmp_path):
        """Test conversion with lineage but no source file specified."""
        adata = self.create_mock_anndata()
        
        project_dir = tmp_path / "test_project_no_source"
        mdv = convert_scanpy_to_mdv(
            folder=str(project_dir),
            scanpy_object=adata,
            track_lineage=True,
            source_file=None  # No source file
        )
        
        # Lineage should still be created
        lineage_file = project_dir / "lineage.json"
        assert lineage_file.exists()
        
        with open(lineage_file) as f:
            lineage = json.load(f)
        
        # But no source files recorded
        assert len(lineage['source_files']) == 0
    
    def test_get_lineage_from_mdvproject(self, tmp_path):
        """Test querying lineage from MDVProject object."""
        adata = self.create_mock_anndata()
        
        project_dir = tmp_path / "test_project_query"
        mdv = convert_scanpy_to_mdv(
            folder=str(project_dir),
            scanpy_object=adata,
            max_dims=3,
            track_lineage=True
        )
        
        # Query lineage using MDVProject method
        lineage = mdv.get_lineage()
        
        assert lineage is not None
        assert 'umts_version' in lineage
        assert lineage['conversion']['parameters']['max_dims'] == 3
    
    def test_get_lineage_when_not_tracked(self, tmp_path):
        """Test querying lineage when it wasn't tracked."""
        adata = self.create_mock_anndata()
        
        project_dir = tmp_path / "test_project_no_lineage"
        mdv = convert_scanpy_to_mdv(
            folder=str(project_dir),
            scanpy_object=adata,
            track_lineage=False
        )
        
        # Should return None
        lineage = mdv.get_lineage()
        assert lineage is None
    
    def test_backward_compatibility(self, tmp_path):
        """Test that old code without lineage params still works."""
        adata = self.create_mock_anndata()
        
        project_dir = tmp_path / "test_project_compat"
        
        # Call without new parameters (should work with defaults)
        mdv = convert_scanpy_to_mdv(
            folder=str(project_dir),
            scanpy_object=adata,
            max_dims=3,
            delete_existing=True
        )
        
        # Should still create project successfully
        assert os.path.exists(str(project_dir / "datafile.h5"))
        assert os.path.exists(str(project_dir / "datasources.json"))
        
        # Lineage should be created by default (track_lineage=True)
        lineage_file = project_dir / "lineage.json"
        assert lineage_file.exists()
    
    def test_lineage_with_label_parameter(self, tmp_path):
        """Test lineage tracking with label parameter."""
        adata = self.create_mock_anndata()
        
        project_dir = tmp_path / "test_project_label"
        mdv = convert_scanpy_to_mdv(
            folder=str(project_dir),
            scanpy_object=adata,
            label="rna_",
            track_lineage=True
        )
        
        lineage = mdv.get_lineage()
        assert lineage['conversion']['parameters']['label'] == "rna_"
    
    def test_reopen_project_and_query_lineage(self, tmp_path):
        """Test reopening a project and querying its lineage."""
        adata = self.create_mock_anndata()
        
        project_dir = tmp_path / "test_project_reopen"
        
        # Create project with lineage
        mdv1 = convert_scanpy_to_mdv(
            folder=str(project_dir),
            scanpy_object=adata,
            max_dims=2,
            track_lineage=True
        )
        
        original_lineage = mdv1.get_lineage()
        
        # Close and reopen project
        mdv2 = MDVProject(str(project_dir))
        
        # Should be able to query lineage
        reopened_lineage = mdv2.get_lineage()
        
        assert reopened_lineage is not None
        assert reopened_lineage['conversion']['parameters']['max_dims'] == 2
        assert reopened_lineage['umts_version'] == original_lineage['umts_version']

