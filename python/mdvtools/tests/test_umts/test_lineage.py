"""
Unit tests for UMTS lineage tracking functionality.
"""

import pytest
import os
import json
import tempfile
import hashlib
from pathlib import Path

from mdvtools.umts import (
    LineageTracker,
    compute_file_hash,
    get_package_versions,
    get_timestamp,
    get_file_metadata,
    get_environment_info
)


class TestUtils:
    """Test utility functions."""
    
    def test_compute_file_hash(self, tmp_path):
        """Test file hash computation."""
        # Create a test file
        test_file = tmp_path / "test.txt"
        test_content = b"Hello, UMTS!"
        test_file.write_bytes(test_content)
        
        # Compute expected hash
        expected_hash = hashlib.sha256(test_content).hexdigest()
        
        # Test hash computation
        computed_hash = compute_file_hash(str(test_file))
        assert computed_hash == expected_hash
    
    def test_compute_file_hash_nonexistent(self):
        """Test hash computation with nonexistent file."""
        with pytest.raises(FileNotFoundError):
            compute_file_hash("/nonexistent/file.txt")
    
    def test_get_package_versions(self):
        """Test package version detection."""
        versions = get_package_versions(['numpy', 'pandas'])
        
        assert isinstance(versions, dict)
        assert 'numpy' in versions
        assert 'pandas' in versions
        
        # Check that numpy is installed (should be)
        assert versions['numpy'] != "not installed"
    
    def test_get_package_versions_default(self):
        """Test default package list."""
        versions = get_package_versions()
        
        # Should include default packages
        assert 'numpy' in versions
        assert 'pandas' in versions
        assert 'scipy' in versions
    
    def test_get_timestamp(self):
        """Test timestamp generation."""
        timestamp = get_timestamp()
        
        # Should be in ISO8601 format
        assert isinstance(timestamp, str)
        assert 'T' in timestamp
        assert timestamp.endswith('Z')
        
        # Should be parseable
        from datetime import datetime
        dt = datetime.fromisoformat(timestamp.replace('Z', '+00:00'))
        assert dt is not None
    
    def test_get_file_metadata(self, tmp_path):
        """Test file metadata extraction."""
        # Create a test file
        test_file = tmp_path / "metadata_test.txt"
        test_content = b"Test content"
        test_file.write_bytes(test_content)
        
        metadata = get_file_metadata(str(test_file))
        
        assert 'path' in metadata
        assert 'size_bytes' in metadata
        assert 'modified_timestamp' in metadata
        
        assert metadata['size_bytes'] == len(test_content)
        assert metadata['path'] == str(test_file.absolute())
    
    def test_get_environment_info(self):
        """Test environment information gathering."""
        env_info = get_environment_info()
        
        assert 'python_version' in env_info
        assert 'platform' in env_info
        assert 'packages' in env_info
        
        assert isinstance(env_info['packages'], dict)


class TestLineageTracker:
    """Test LineageTracker class."""
    
    def test_initialization(self):
        """Test LineageTracker initialization."""
        tracker = LineageTracker()
        
        data = tracker.to_dict()
        assert 'umts_version' in data
        assert 'created_timestamp' in data
        assert 'source_files' in data
        assert 'conversion' in data
        assert 'environment' in data
        assert 'notes' in data
        
        assert isinstance(data['source_files'], list)
        assert len(data['source_files']) == 0
    
    def test_record_source(self, tmp_path):
        """Test recording source file."""
        tracker = LineageTracker()
        
        # Create a test file
        test_file = tmp_path / "source.h5ad"
        test_file.write_bytes(b"Mock AnnData")
        
        tracker.record_source(str(test_file))
        
        data = tracker.to_dict()
        assert len(data['source_files']) == 1
        
        source = data['source_files'][0]
        assert 'path' in source
        assert 'sha256' in source
        assert 'size_bytes' in source
        assert source['size_bytes'] == 12  # len(b"Mock AnnData")
    
    def test_record_source_no_hash(self, tmp_path):
        """Test recording source without computing hash."""
        tracker = LineageTracker()
        
        test_file = tmp_path / "source.h5ad"
        test_file.write_bytes(b"Mock AnnData")
        
        tracker.record_source(str(test_file), compute_hash=False)
        
        data = tracker.to_dict()
        source = data['source_files'][0]
        assert source['sha256'] == "not computed"
    
    def test_record_source_nonexistent(self):
        """Test recording nonexistent source file."""
        tracker = LineageTracker()
        
        with pytest.raises(FileNotFoundError):
            tracker.record_source("/nonexistent/file.h5ad")
    
    def test_record_parameters(self):
        """Test recording conversion parameters."""
        tracker = LineageTracker()
        
        params = {
            'max_dims': 3,
            'delete_existing': True,
            'chunk_data': False
        }
        
        tracker.record_parameters(params, function_name='convert_scanpy_to_mdv')
        
        data = tracker.to_dict()
        assert data['conversion']['function'] == 'convert_scanpy_to_mdv'
        assert data['conversion']['parameters'] == params
    
    def test_record_environment(self):
        """Test recording environment information."""
        tracker = LineageTracker()
        tracker.record_environment()
        
        data = tracker.to_dict()
        assert 'python_version' in data['environment']
        assert 'platform' in data['environment']
        assert 'packages' in data['environment']
    
    def test_add_note(self):
        """Test adding notes."""
        tracker = LineageTracker()
        
        tracker.add_note("First note")
        tracker.add_note("Second note")
        
        data = tracker.to_dict()
        assert len(data['notes']) == 2
        
        assert data['notes'][0]['message'] == "First note"
        assert data['notes'][1]['message'] == "Second note"
        
        # Check timestamps are present
        assert 'timestamp' in data['notes'][0]
        assert 'timestamp' in data['notes'][1]
    
    def test_save_and_load(self, tmp_path):
        """Test saving and loading lineage."""
        tracker = LineageTracker()
        
        # Add some data
        tracker.record_parameters({'max_dims': 3})
        tracker.record_environment()
        tracker.add_note("Test conversion")
        
        # Save
        saved_path = tracker.save(str(tmp_path))
        
        assert os.path.exists(saved_path)
        assert saved_path == str(tmp_path / "lineage.json")
        
        # Load
        loaded_data = LineageTracker.load(str(tmp_path))
        
        assert loaded_data is not None
        assert loaded_data['conversion']['parameters']['max_dims'] == 3
        assert len(loaded_data['notes']) == 1
        assert loaded_data['notes'][0]['message'] == "Test conversion"
    
    def test_load_nonexistent(self, tmp_path):
        """Test loading from directory without lineage file."""
        loaded = LineageTracker.load(str(tmp_path))
        assert loaded is None
    
    def test_from_dict(self):
        """Test creating LineageTracker from dictionary."""
        data = {
            'umts_version': '0.1.0',
            'created_timestamp': '2025-01-01T00:00:00Z',
            'source_files': [],
            'conversion': {'function': 'test', 'parameters': {}},
            'environment': {},
            'notes': []
        }
        
        tracker = LineageTracker.from_dict(data)
        
        restored_data = tracker.to_dict()
        assert restored_data['umts_version'] == '0.1.0'
        assert restored_data['conversion']['function'] == 'test'
    
    def test_save_creates_directory(self, tmp_path):
        """Test that save creates directory if it doesn't exist."""
        tracker = LineageTracker()
        
        new_dir = tmp_path / "new_project"
        tracker.save(str(new_dir))
        
        assert new_dir.exists()
        assert (new_dir / "lineage.json").exists()


class TestLineageIntegration:
    """Integration tests for lineage tracking."""
    
    def test_full_workflow(self, tmp_path):
        """Test complete lineage tracking workflow."""
        # Create mock source file
        source_file = tmp_path / "input.h5ad"
        source_file.write_bytes(b"Mock single-cell data")
        
        # Initialize tracker
        tracker = LineageTracker()
        
        # Record source
        tracker.record_source(str(source_file))
        
        # Record conversion parameters
        tracker.record_parameters({
            'max_dims': 3,
            'delete_existing': True,
            'label': '',
            'chunk_data': False
        }, function_name='convert_scanpy_to_mdv')
        
        # Record environment
        tracker.record_environment()
        
        # Add note
        tracker.add_note("Test conversion run")
        
        # Save to project directory
        project_dir = tmp_path / "test_project"
        tracker.save(str(project_dir))
        
        # Verify file exists and is valid JSON
        lineage_file = project_dir / "lineage.json"
        assert lineage_file.exists()
        
        with open(lineage_file) as f:
            data = json.load(f)
        
        # Verify all expected fields are present
        assert data['umts_version'] == '0.1.0'
        assert len(data['source_files']) == 1
        assert data['conversion']['function'] == 'convert_scanpy_to_mdv'
        assert 'python_version' in data['environment']
        assert len(data['notes']) == 1

