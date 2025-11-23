"""
Lineage tracking for UMTS (Universal Modality Translator System).

Provides the LineageTracker class for recording and querying data provenance
information throughout the conversion pipeline.
"""

import json
import os
from typing import Dict, List, Optional, Any
from .utils import (
    compute_file_hash,
    get_file_metadata,
    get_environment_info,
    get_timestamp
)


# UMTS version - increment when making breaking changes to lineage format
UMTS_VERSION = "0.1.0"


class LineageTracker:
    """
    Tracks lineage and provenance information for data conversions.
    
    Records:
    - Source files with hashes and metadata
    - Conversion function and parameters
    - Environment information (Python version, package versions)
    - User annotations and notes
    - Timestamps
    
    Example:
        tracker = LineageTracker()
        tracker.record_source("input.h5ad")
        tracker.record_parameters({
            "max_dims": 3,
            "delete_existing": True
        })
        tracker.record_environment()
        tracker.add_note("Initial conversion from Seurat export")
        tracker.save("/path/to/mdv/project")
    """
    
    def __init__(self):
        """Initialize a new LineageTracker."""
        self.data = {
            'umts_version': UMTS_VERSION,
            'created_timestamp': get_timestamp(),
            'source_files': [],
            'conversion': {
                'function': None,
                'parameters': {}
            },
            'environment': {},
            'notes': []
        }
    
    def record_source(self, file_path: str, compute_hash: bool = True) -> None:
        """
        Record a source file with its metadata.
        
        Args:
            file_path: Path to the source file
            compute_hash: Whether to compute SHA256 hash (default: True).
                         Set to False for very large files to save time.
                         
        Raises:
            FileNotFoundError: If the file doesn't exist
        """
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"Source file not found: {file_path}")
        
        file_info = get_file_metadata(file_path)
        
        if compute_hash:
            try:
                file_info['sha256'] = compute_file_hash(file_path)
            except Exception as e:
                file_info['sha256'] = f"error: {str(e)}"
                file_info['sha256_error'] = True
        else:
            file_info['sha256'] = "not computed"
        
        self.data['source_files'].append(file_info)
    
    def record_parameters(
        self,
        params_dict: Dict[str, Any],
        function_name: Optional[str] = None
    ) -> None:
        """
        Record conversion function parameters.
        
        Args:
            params_dict: Dictionary of parameter names and values
            function_name: Name of the conversion function (optional)
        """
        self.data['conversion']['parameters'] = params_dict.copy()
        
        if function_name:
            self.data['conversion']['function'] = function_name
    
    def record_environment(self) -> None:
        """
        Record information about the current Python environment.
        
        Captures Python version, platform, and versions of key packages.
        """
        self.data['environment'] = get_environment_info()
    
    def add_note(self, message: str) -> None:
        """
        Add a user annotation or note to the lineage.
        
        Args:
            message: The note text to add
        """
        self.data['notes'].append({
            'timestamp': get_timestamp(),
            'message': message
        })
    
    def to_dict(self) -> Dict:
        """
        Convert lineage data to dictionary format.
        
        Returns:
            Dictionary containing all lineage information
        """
        return self.data.copy()
    
    def save(self, project_dir: str, filename: str = "lineage.json") -> str:
        """
        Save lineage information to a JSON file.
        
        Args:
            project_dir: Directory to save the lineage file in
            filename: Name of the lineage file (default: "lineage.json")
            
        Returns:
            Path to the saved lineage file
            
        Raises:
            IOError: If there's an error writing the file
        """
        if not os.path.exists(project_dir):
            os.makedirs(project_dir, exist_ok=True)
        
        lineage_path = os.path.join(project_dir, filename)
        
        try:
            with open(lineage_path, 'w') as f:
                json.dump(self.data, f, indent=2)
        except Exception as e:
            raise IOError(f"Error saving lineage file: {str(e)}")
        
        return lineage_path
    
    @staticmethod
    def load(project_dir: str, filename: str = "lineage.json") -> Optional[Dict]:
        """
        Load lineage information from a JSON file.
        
        Args:
            project_dir: Directory containing the lineage file
            filename: Name of the lineage file (default: "lineage.json")
            
        Returns:
            Dictionary containing lineage information, or None if file doesn't exist
            
        Raises:
            json.JSONDecodeError: If the file contains invalid JSON
            IOError: If there's an error reading the file
        """
        lineage_path = os.path.join(project_dir, filename)
        
        if not os.path.exists(lineage_path):
            return None
        
        try:
            with open(lineage_path, 'r') as f:
                return json.load(f)
        except json.JSONDecodeError as e:
            raise json.JSONDecodeError(
                f"Invalid JSON in lineage file: {str(e)}",
                e.doc,
                e.pos
            )
        except Exception as e:
            raise IOError(f"Error reading lineage file: {str(e)}")
    
    @staticmethod
    def from_dict(data: Dict) -> 'LineageTracker':
        """
        Create a LineageTracker instance from a dictionary.
        
        Args:
            data: Dictionary containing lineage data
            
        Returns:
            New LineageTracker instance with the provided data
        """
        tracker = LineageTracker()
        tracker.data = data.copy()
        return tracker

