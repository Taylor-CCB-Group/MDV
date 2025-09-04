"""
Flask routes for metadata extraction functionality.

This module contains all the route handlers for fetching metadata from
various dataset formats including Zarr stores and SpatialData datasets.
"""

import os
import asyncio
from flask import request, jsonify
from mdvtools.zarr_utils import get_zarr_extractor
from mdvtools.logging_config import get_logger

logger = get_logger(__name__)

def register_metadata_routes(project_bp, project):
    """
    Register metadata-related routes with the project blueprint.
    
    Args:
        project_bp: The Flask blueprint to register routes with
        project: The MDVProject instance
    """
    
    @project_bp.route("/get_metadata", methods=["GET"])
    def get_metadata():
        """
        Fetch metadata from a dataset URL (supports Zarr and SpatialData formats).
        
        Query Parameters:
            url: The URL of the dataset to fetch metadata from
            
        Returns:
            JSON response with metadata structure or error message
        """
        logger.info("=== METADATA REQUEST ===")
        logger.info(f"Request method: {request.method}")
        logger.info(f"Request path: {request.path}")
        logger.info(f"Request args: {request.args}")
        logger.info(f"Project: {project.id}")
        
        dataset_url = None
        try:
            dataset_url = request.args.get('url')
            logger.info(f"Extracted dataset_url: {dataset_url}")
            
            if not dataset_url:
                logger.info("ERROR: Missing 'url' parameter")
                return jsonify({
                    "error": "Missing 'url' parameter",
                    "message": "Please provide a 'url' query parameter with the dataset URL"
                }), 400
                
            # Validate URL format
            if not dataset_url.startswith(('http://', 'https://')):
                logger.info(f"ERROR: Invalid URL format: {dataset_url}")
                return jsonify({
                    "error": "Invalid URL format", 
                    "message": "URL must start with http:// or https://"
                }), 400
            
            # Convert localhost to host.docker.internal if running in Docker
            if 'localhost' in dataset_url and os.path.exists('/.dockerenv'):
                dataset_url = dataset_url.replace('localhost', 'host.docker.internal')
                logger.info(f"Running in Docker, converted URL to: {dataset_url}")
            
            logger.info(f"URL validation passed, fetching metadata for: {dataset_url}")
            
            # Extract metadata using zarr utilities (supports both Zarr and SpatialData)
            logger.info("Getting metadata extractor...")
            extractor = get_zarr_extractor()
            logger.info(f"Metadata extractor: {extractor}")
            
            logger.info("Calling fetch_zarr_metadata...")
            print("Creating new event loop for async metadata call...")
            # Run the async function in the current event loop
            loop = asyncio.new_event_loop()
            asyncio.set_event_loop(loop)
            try:
                print("Running fetch_zarr_metadata in event loop...")
                metadata = loop.run_until_complete(extractor.fetch_zarr_metadata(dataset_url))
                print(f"Event loop completed, metadata type: {type(metadata)}")
            except Exception as loop_e:
                print(f"Exception in event loop: {loop_e}")
                raise loop_e
            finally:
                print("Closing event loop...")
                loop.close()
                print("Event loop closed")
            logger.info(f"Metadata extraction successful, keys: {list(metadata.keys()) if isinstance(metadata, dict) else type(metadata)}")
            print(f"Metadata extraction successful, keys: {list(metadata.keys()) if isinstance(metadata, dict) else type(metadata)}")
            
            response_data = {
                "success": True,
                "metadata": metadata,
                "url": dataset_url
            }
            logger.info("Returning successful response")
            return jsonify(response_data)
            
        except Exception as e:
            logger.info(f"EXCEPTION in get_metadata: {str(e)}")
            logger.info(f"Exception type: {type(e)}")
            import traceback
            logger.info(f"Full traceback: {traceback.format_exc()}")
            
            return jsonify({
                "error": "Metadata extraction failed",
                "message": str(e),
                "url": dataset_url
            }), 500