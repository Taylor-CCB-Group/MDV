"""
Metadata extraction and processing module for MDV.

This module provides functionality for fetching and processing metadata
from various dataset formats including Zarr stores and SpatialData datasets.
"""

from .routes import register_metadata_routes

__all__ = ['register_metadata_routes']