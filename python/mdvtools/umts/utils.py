"""
Utility functions for UMTS (Universal Modality Translator System).

Provides helper functions for file hashing, package version detection,
and timestamp generation for lineage tracking.
"""

import hashlib
import os
import platform
import sys
from datetime import datetime, timezone
from typing import Dict, Optional
import importlib.metadata


def compute_file_hash(filepath: str, algorithm: str = "sha256") -> str:
    """
    Compute cryptographic hash of a file.
    
    Args:
        filepath: Path to the file to hash
        algorithm: Hash algorithm to use (default: sha256)
        
    Returns:
        Hexadecimal string representation of the hash
        
    Raises:
        FileNotFoundError: If the file doesn't exist
        IOError: If there's an error reading the file
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")
    
    hash_obj = hashlib.new(algorithm)
    
    # Read file in chunks to handle large files efficiently
    chunk_size = 8192
    with open(filepath, 'rb') as f:
        while chunk := f.read(chunk_size):
            hash_obj.update(chunk)
    
    return hash_obj.hexdigest()


def get_package_versions(packages: Optional[list] = None) -> Dict[str, str]:
    """
    Get versions of specified packages.
    
    Args:
        packages: List of package names to check. If None, checks default
                 packages relevant to MDV/UMTS (mdvtools, scanpy, scipy, 
                 numpy, pandas, anndata, h5py).
                 
    Returns:
        Dictionary mapping package names to version strings.
        If a package is not found, its value will be "not installed".
    """
    if packages is None:
        packages = [
            'mdvtools',
            'scanpy', 
            'scipy',
            'numpy',
            'pandas',
            'anndata',
            'h5py',
            'polars'
        ]
    
    versions = {}
    for package in packages:
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = "not installed"
    
    return versions


def get_timestamp() -> str:
    """
    Get current timestamp in ISO8601 format with UTC timezone.
    
    Returns:
        ISO8601 formatted timestamp string (e.g., "2025-11-23T10:30:00Z")
    """
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def get_file_metadata(filepath: str) -> Dict[str, any]:
    """
    Get metadata about a file.
    
    Args:
        filepath: Path to the file
        
    Returns:
        Dictionary containing:
            - path: absolute file path
            - size_bytes: file size in bytes
            - modified_timestamp: last modification time in ISO8601 format
            
    Raises:
        FileNotFoundError: If the file doesn't exist
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")
    
    stat = os.stat(filepath)
    modified_dt = datetime.fromtimestamp(stat.st_mtime, tz=timezone.utc)
    
    return {
        'path': os.path.abspath(filepath),
        'size_bytes': stat.st_size,
        'modified_timestamp': modified_dt.strftime("%Y-%m-%dT%H:%M:%SZ")
    }


def get_environment_info() -> Dict[str, any]:
    """
    Get information about the current Python environment.
    
    Returns:
        Dictionary containing:
            - python_version: Python version string
            - platform: Platform/OS information
            - packages: Dictionary of package versions
    """
    return {
        'python_version': sys.version.split()[0],
        'platform': platform.platform(),
        'packages': get_package_versions()
    }

