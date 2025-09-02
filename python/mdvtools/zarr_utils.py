"""
Zarr utilities for MDV platform - backend implementation of Zarr metadata extraction
"""

import zarr
import aiohttp
import asyncio
import json
from typing import Dict, Any, Optional, List
from urllib.parse import urlparse, urljoin
import logging

logger = logging.getLogger(__name__)


class ZarrMetadataExtractor:
    """Extracts metadata from Zarr stores and SpatialData formats for the MDV platform"""

    def __init__(self):
        self.timeout = aiohttp.ClientTimeout(total=30)

    async def fetch_zarr_metadata(self, url: str) -> Dict[str, Any]:
        """
        Fetch metadata from a Zarr dataset or SpatialData URL using async HTTP requests.
        Automatically detects the format and extracts appropriate metadata.
        
        Args:
            url: URL to the Zarr dataset (can be HTTP/HTTPS)
            
        Returns:
            Dictionary containing structured metadata for frontend display
            
        Raises:
            Exception: If metadata extraction fails
        """
        try:
            print(f"[ZARR] Starting fetch_zarr_metadata for: {url}")
            logger.info(f"Fetching Zarr metadata from: {url}")
            
            if not url.startswith(('http://', 'https://')):
                raise ValueError(f"Unsupported URL scheme: {url}")
            
            print(f"[ZARR] Creating aiohttp session...")
            # Use async HTTP approach to avoid fsspec conflicts
            async with aiohttp.ClientSession(timeout=self.timeout) as session:
                print(f"[ZARR] Session created, detecting dataset format...")
                dataset_format = await self._detect_dataset_format(session, url)
                print(f"[ZARR] Detected format: {dataset_format}")
                
                if dataset_format == "spatialdata":
                    print(f"[ZARR] Processing as SpatialData format...")
                    metadata = await self._fetch_spatialdata_metadata(session, url)
                else:
                    print(f"[ZARR] Processing as standard Zarr format...")
                    metadata = await self._fetch_http_zarr_metadata(session, url)
                    
                print(f"[ZARR] Metadata extraction returned: {type(metadata)}")
            
            print(f"[ZARR] Successfully extracted metadata")
            logger.info("Successfully extracted metadata")
            return metadata
            
        except Exception as e:
            logger.error(f"Error fetching Zarr metadata from {url}: {str(e)}")
            raise Exception(f"Failed to fetch metadata from {url}: {str(e)}")

    async def _detect_dataset_format(self, session: aiohttp.ClientSession, url: str) -> str:
        """
        Detect if the dataset is SpatialData format or standard Zarr format.
        SpatialData typically has organized groups: images, points, shapes, tables.
        """
        if not url.endswith('/'):
            url = url + '/'
        
        logger.info(f"Detecting dataset format for: {url}")
        
        # Check for SpatialData indicators
        spatialdata_groups = ['images', 'points', 'shapes', 'tables']
        found_groups = 0
        
        for group in spatialdata_groups:
            try:
                group_url = urljoin(url, f'{group}/')
                zgroup_url = urljoin(group_url, '.zgroup')
                
                async with session.get(zgroup_url) as response:
                    if response.status == 200:
                        found_groups += 1
                        logger.info(f"Found SpatialData group: {group}")
            except Exception as e:
                logger.debug(f"Could not check group {group}: {e}")
        
        # If we found 2 or more SpatialData-specific groups, classify as SpatialData
        if found_groups >= 2:
            logger.info("Detected SpatialData format")
            return "spatialdata"
        else:
            logger.info("Detected standard Zarr format")
            return "zarr"

    async def _fetch_spatialdata_metadata(self, session: aiohttp.ClientSession, url: str) -> Dict[str, Any]:
        """
        Fetch metadata from a SpatialData-formatted Zarr store.
        """
        logger.info(f"Fetching SpatialData metadata from: {url}")
        
        if not url.endswith('/'):
            url = url + '/'
        
        metadata = {
            "datasetFormat": "spatialdata",
            "datasetStructure": {
                "groups": [],
                "arrays": []
            },
            "spatialData": {
                "coordinateSystems": {},
                "elements": {},
                "transformations": []
            },
            "images": {},
            "tables": {},
            "labels": {},
            "points": {},
            "shapes": {},
            "groupAttributes": {},
            "rawStructure": True
        }
        
        # Get root attributes
        try:
            attrs_url = urljoin(url, '.zattrs')
            async with session.get(attrs_url) as response:
                if response.status == 200:
                    attrs_text = await response.text()
                    metadata["groupAttributes"] = json.loads(attrs_text)
                    logger.info(f"Found root attributes: {list(metadata['groupAttributes'].keys())}")
        except Exception as e:
            logger.warning(f"Could not fetch root attributes: {e}")
        
        # Process SpatialData groups
        spatialdata_groups = {
            'images': self._process_images_group,
            'points': self._process_points_group, 
            'shapes': self._process_shapes_group,
            'tables': self._process_tables_group
        }
        
        for group_name, processor in spatialdata_groups.items():
            try:
                group_url = urljoin(url, f'{group_name}/')
                zgroup_url = urljoin(group_url, '.zgroup')
                
                async with session.get(zgroup_url) as response:
                    if response.status == 200:
                        logger.info(f"Processing SpatialData group: {group_name}")
                        metadata["datasetStructure"]["groups"].append(group_name)
                        await processor(session, group_url, metadata)
                        
            except Exception as e:
                logger.info(f"Group {group_name} not found or inaccessible: {e}")
        
        logger.info(f"SpatialData metadata extraction completed. Found {len(metadata['datasetStructure']['groups'])} groups, {len(metadata['datasetStructure']['arrays'])} arrays")
        return metadata

    async def _fetch_http_zarr_metadata(self, session: aiohttp.ClientSession, url: str) -> Dict[str, Any]:
        """
        Fetch Zarr metadata over HTTP using async requests.
        """
        print(f"[ZARR] _fetch_http_zarr_metadata called with url: {url}")
        logger.info(f"Fetching HTTP Zarr metadata from: {url}")
        
        # Ensure URL ends with /
        if not url.endswith('/'):
            url = url + '/'
        
        print(f"[ZARR] Normalized URL: {url}")
        
        # Try to fetch consolidated metadata first
        metadata_url = urljoin(url, '.zmetadata')
        
        try:
            print(f"[ZARR] Attempting to fetch consolidated metadata from: {metadata_url}")
            logger.info(f"Attempting to fetch consolidated metadata from: {metadata_url}")
            async with session.get(metadata_url) as response:
                print(f"[ZARR] Consolidated metadata response status: {response.status}")
                if response.status == 200:
                    logger.info("Found consolidated metadata")
                    consolidated_meta = await response.json()
                    print(f"[ZARR] Consolidated metadata keys: {list(consolidated_meta.keys())}")
                    return await self._parse_consolidated_metadata(consolidated_meta, url)
                else:
                    print(f"[ZARR] No consolidated metadata found (status: {response.status}), trying direct exploration")
                    logger.info(f"No consolidated metadata found (status: {response.status}), trying direct exploration")
                    
        except Exception as e:
            print(f"[ZARR] Exception fetching consolidated metadata: {e}")
            logger.info(f"Failed to fetch consolidated metadata: {e}")
        
        # Fallback: try to explore structure directly
        print(f"[ZARR] Falling back to direct structure exploration")
        return await self._explore_zarr_http_structure(session, url)

    async def _parse_consolidated_metadata(self, consolidated_meta: Dict, url: str) -> Dict[str, Any]:
        """Parse consolidated Zarr metadata."""
        logger.info("Parsing consolidated metadata")
        
        metadata = {
            "datasetFormat": "zarr",
            "datasetStructure": {
                "groups": [],
                "arrays": []
            },
            "spatialData": None,
            "images": {},
            "tables": {},
            "labels": {},
            "groupAttributes": {},
            "rawStructure": True
        }
        
        # Extract root attributes from metadata
        if 'metadata' in consolidated_meta and '' in consolidated_meta['metadata']:
            root_meta = consolidated_meta['metadata']['']
            if 'attributes' in root_meta:
                metadata["groupAttributes"] = root_meta['attributes']
        
        # Process all items in the zarr structure
        zarr_meta = consolidated_meta.get('metadata', {})
        
        for path, item_meta in zarr_meta.items():
            if not path:  # Skip root
                continue
                
            if 'zarr_format' in item_meta:  # This is an array
                metadata["datasetStructure"]["arrays"].append(path)
                
                # Classify the array
                array_name = path.split('/')[-1]
                if self._is_image_name(array_name):
                    metadata["images"][array_name] = self._extract_array_metadata_from_meta(item_meta)
                elif self._is_table_name(array_name):
                    metadata["tables"][array_name] = self._extract_array_metadata_from_meta(item_meta)
                elif self._is_label_name(array_name):
                    metadata["labels"][array_name] = self._extract_array_metadata_from_meta(item_meta)
                    
            elif 'attributes' in item_meta:  # This is a group
                metadata["datasetStructure"]["groups"].append(path)
        
        logger.info(f"Parsed consolidated metadata: {len(metadata['datasetStructure']['arrays'])} arrays, {len(metadata['datasetStructure']['groups'])} groups")
        return metadata

    async def _explore_zarr_http_structure(self, session: aiohttp.ClientSession, url: str) -> Dict[str, Any]:
        """Explore Zarr structure via direct HTTP requests."""
        logger.info("Exploring Zarr structure via HTTP")
        
        metadata = {
            "datasetFormat": "zarr",
            "datasetStructure": {
                "groups": [],
                "arrays": []
            },
            "spatialData": None,
            "images": {},
            "tables": {},
            "labels": {},
            "groupAttributes": {},
            "rawStructure": True
        }
        
        # Try to get root attributes
        try:
            attrs_url = urljoin(url, '.zattrs')
            async with session.get(attrs_url) as response:
                if response.status == 200:
                    attrs_text = await response.text()
                    metadata["groupAttributes"] = json.loads(attrs_text)
                    logger.info(f"Found root attributes: {list(metadata['groupAttributes'].keys())}")
        except Exception as e:
            logger.warning(f"Could not fetch root attributes: {e}")
        
        # Look for known Xenium structures
        xenium_structures = ['cells.zarr/', 'transcripts.zarr/', 'analysis.zarr/', 'cells.zarr.zip', 'transcripts.zarr.zip']
        
        for structure in xenium_structures:
            try:
                structure_url = urljoin(url, structure)
                logger.info(f"Checking for structure: {structure_url}")
                
                # Check if it's a zarr group by looking for .zgroup
                zgroup_url = urljoin(structure_url, '.zgroup')
                async with session.get(zgroup_url) as response:
                    if response.status == 200:
                        logger.info(f"Found Zarr group: {structure}")
                        metadata["datasetStructure"]["groups"].append(structure.rstrip('/'))
                        
                        # Get attributes for this group
                        attrs_url = urljoin(structure_url, '.zattrs')
                        try:
                            async with session.get(attrs_url) as response:
                                if response.status == 200:
                                    # Handle JSON parsing regardless of content type
                                    attrs_text = await response.text()
                                    group_attrs = json.loads(attrs_text)
                                    logger.info(f"Group {structure} attributes: {list(group_attrs.keys())}")
                                    
                                    # Store attributes in groupAttributes with path prefix
                                    group_name = structure.rstrip('/').replace('.zarr', '')
                                    metadata["groupAttributes"][group_name] = group_attrs
                        except Exception as e:
                            logger.warning(f"Could not fetch attributes for {structure}: {e}")
                        
                        # Explore arrays within this group
                        await self._explore_group_structure(session, structure_url, structure.rstrip('/'), metadata)
                        
            except Exception as e:
                logger.info(f"Structure {structure} not found or inaccessible: {e}")
                continue
        
        logger.info(f"HTTP-based structure exploration completed. Found {len(metadata['datasetStructure']['groups'])} groups, {len(metadata['datasetStructure']['arrays'])} arrays")
        return metadata

    async def _explore_group_structure(self, session: aiohttp.ClientSession, base_url: str, group_path: str, metadata: Dict[str, Any]):
        """Explore the structure within a Zarr group."""
        logger.info(f"Exploring group structure: {group_path}")
        
        # Known array patterns in Xenium data
        known_arrays = {
            'cells.zarr': ['cell_id', 'cell_summary'],
            'transcripts.zarr': ['codeword_category', 'gene_category', 'metrics_density'],
            'analysis.zarr': ['cell_groups']
        }
        
        group_name = group_path.replace('.zarr', '')
        arrays_to_check = known_arrays.get(f'{group_name}.zarr', [])
        
        for array_name in arrays_to_check:
            try:
                array_url = urljoin(base_url, f'{array_name}/')
                zarray_url = urljoin(array_url, '.zarray')
                
                async with session.get(zarray_url) as response:
                    if response.status == 200:
                        # Handle JSON parsing regardless of content type
                        zarray_text = await response.text()
                        array_meta = json.loads(zarray_text)
                        full_array_path = f"{group_path}/{array_name}"
                        metadata["datasetStructure"]["arrays"].append(full_array_path)
                        
                        logger.info(f"Found array: {full_array_path} with shape {array_meta.get('shape', 'unknown')}")
                        
                        # Classify and store array metadata
                        if self._is_table_name(array_name) or group_name in ['cells', 'analysis']:
                            metadata["tables"][array_name] = {
                                "shape": array_meta.get('shape', []),
                                "dtype": array_meta.get('dtype', 'unknown'),
                                "chunks": array_meta.get('chunks', []),
                                "group": group_name
                            }
                            
                        # Get array attributes
                        try:
                            attrs_url = urljoin(array_url, '.zattrs')
                            async with session.get(attrs_url) as response:
                                if response.status == 200:
                                    attrs_text = await response.text()
                                    array_attrs = json.loads(attrs_text)
                                    metadata["tables"][array_name]["attributes"] = array_attrs
                        except Exception as e:
                            logger.info(f"Could not get attributes for array {array_name}: {e}")
                            
            except Exception as e:
                logger.info(f"Array {array_name} in {group_path} not found: {e}")
                continue

    def _is_image_name(self, name: str) -> bool:
        """Determine if name suggests image data"""
        image_patterns = ['image', 'img', 'tissue', 'dapi', 'he', 'rgb']
        name_lower = name.lower()
        return any(pattern in name_lower for pattern in image_patterns)

    def _is_table_name(self, name: str) -> bool:
        """Determine if name suggests tabular data"""
        table_patterns = ['table', 'obs', 'var', 'cell', 'gene', 'transcript']
        name_lower = name.lower()
        return any(pattern in name_lower for pattern in table_patterns)

    def _is_label_name(self, name: str) -> bool:
        """Determine if name suggests label data"""
        label_patterns = ['label', 'segmentation', 'mask', 'nuclei', 'cell_id']
        name_lower = name.lower()
        return any(pattern in name_lower for pattern in label_patterns)

    def _extract_array_metadata_from_meta(self, item_meta: Dict) -> Dict[str, Any]:
        """Extract array metadata from consolidated metadata."""
        attrs = item_meta.get('attributes', {})
        
        return {
            "shape": item_meta.get('shape', []),
            "dtype": item_meta.get('dtype', 'unknown'),
            "chunks": item_meta.get('chunks', []),
            "attributes": attrs
        }

    def _extract_metadata_structure(self, root, url: str) -> Dict[str, Any]:
        """Extract structured metadata from zarr root group"""
        
        metadata = {
            "datasetStructure": {
                "groups": [],
                "arrays": []
            },
            "spatialData": None,
            "images": {},
            "tables": {},
            "labels": {},
            "groupAttributes": {},
            "rawStructure": True
        }
        
        # Get root attributes
        if hasattr(root, 'attrs'):
            metadata["groupAttributes"] = dict(root.attrs)
        
        # Recursively explore structure
        self._explore_group(root, metadata, "")
        
        # Check for spatial data conventions
        if self._is_spatial_data_format(root):
            metadata["spatialData"] = self._extract_spatial_data_info(root)
        
        return metadata

    def _explore_group(self, group, metadata: Dict[str, Any], path_prefix: str):
        """Recursively explore zarr group structure"""
        
        for key in group.keys():
            full_path = f"{path_prefix}/{key}" if path_prefix else key
            
            try:
                item = group[key]
                
                if isinstance(item, zarr.Group):
                    metadata["datasetStructure"]["groups"].append(full_path)
                    # Recurse into subgroups
                    self._explore_group(item, metadata, full_path)
                    
                elif isinstance(item, zarr.Array):
                    metadata["datasetStructure"]["arrays"].append(full_path)
                    
                    # Classify arrays based on name patterns and attributes
                    if self._is_image_array(key, item):
                        metadata["images"][key] = self._extract_image_metadata(item)
                    elif self._is_table_array(key, item):
                        metadata["tables"][key] = self._extract_table_metadata(item)
                    elif self._is_label_array(key, item):
                        metadata["labels"][key] = self._extract_label_metadata(item)
                        
            except Exception as e:
                logger.warning(f"Could not process item {full_path}: {e}")
                continue

    def _is_spatial_data_format(self, root) -> bool:
        """Check if this follows SpatialData conventions"""
        attrs = getattr(root, 'attrs', {})
        return ('spatialdata' in str(attrs).lower() or 
                any(key in attrs for key in ['coordinate_systems', 'elements']))

    def _extract_spatial_data_info(self, root) -> Dict[str, Any]:
        """Extract spatial data specific information"""
        attrs = dict(getattr(root, 'attrs', {}))
        
        spatial_info = {
            "coordinateSystems": {},
            "elements": {},
            "transformations": []
        }
        
        # Look for coordinate system info
        if 'coordinate_systems' in attrs:
            spatial_info["coordinateSystems"] = attrs['coordinate_systems']
        
        # Look for element mappings
        if 'elements' in attrs:
            spatial_info["elements"] = attrs['elements']
            
        return spatial_info

    def _is_image_array(self, name: str, array: zarr.Array) -> bool:
        """Determine if array represents image data"""
        # Check shape - images typically have 2+ dimensions
        if len(array.shape) < 2:
            return False
        
        # Check name patterns
        image_patterns = ['image', 'img', 'tissue', 'dapi', 'he', 'rgb']
        name_lower = name.lower()
        
        return any(pattern in name_lower for pattern in image_patterns)

    def _is_table_array(self, name: str, array: zarr.Array) -> bool:
        """Determine if array represents tabular data"""
        table_patterns = ['table', 'obs', 'var', 'cell', 'gene', 'transcript']
        name_lower = name.lower()
        
        return any(pattern in name_lower for pattern in table_patterns)

    def _is_label_array(self, name: str, array: zarr.Array) -> bool:
        """Determine if array represents label/segmentation data"""
        label_patterns = ['label', 'segmentation', 'mask', 'nuclei', 'cell_id']
        name_lower = name.lower()
        
        return any(pattern in name_lower for pattern in label_patterns)

    def _extract_image_metadata(self, array: zarr.Array) -> Dict[str, Any]:
        """Extract metadata specific to image arrays"""
        attrs = dict(getattr(array, 'attrs', {}))
        
        return {
            "shape": list(array.shape),
            "dtype": str(array.dtype),
            "chunks": list(array.chunks) if hasattr(array, 'chunks') else [],
            "dimensions": attrs.get('dimensions', [f"dim_{i}" for i in range(len(array.shape))]),
            "physicalSizeX": attrs.get('physicalSizeX'),
            "physicalSizeY": attrs.get('physicalSizeY'),
            "physicalSizeZ": attrs.get('physicalSizeZ'),
            "units": attrs.get('units'),
            "channels": attrs.get('channels', [])
        }

    def _extract_table_metadata(self, array: zarr.Array) -> Dict[str, Any]:
        """Extract metadata specific to table arrays"""
        attrs = dict(getattr(array, 'attrs', {}))
        
        # Try to get a small sample of data for preview
        obs_names = []
        var_names = []
        
        try:
            # For 2D arrays, try to interpret as observation x variable matrix
            if len(array.shape) == 2:
                obs_names = [f"obs_{i}" for i in range(min(3, array.shape[0]))]
                var_names = [f"var_{i}" for i in range(min(3, array.shape[1]))]
        except Exception:
            pass
        
        return {
            "shape": list(array.shape),
            "obsNames": obs_names,
            "varNames": var_names,
            "observations": attrs.get('observations', {}),
            "variables": attrs.get('variables', {}),
            "embeddings": attrs.get('embeddings', [])
        }

    def _extract_label_metadata(self, array: zarr.Array) -> Dict[str, Any]:
        """Extract metadata specific to label arrays"""
        attrs = dict(getattr(array, 'attrs', {}))
        
        # Try to get unique values (sample)
        label_values = []
        try:
            if array.size < 1000000:  # Only sample if reasonably sized
                sample_data = array[...]
                unique_vals = list(set(sample_data.flat))[:10]  # Get first 10 unique values
                label_values = [int(v) for v in unique_vals if isinstance(v, (int, float))]
        except Exception:
            pass
        
        return {
            "shape": list(array.shape),
            "dtype": str(array.dtype),
            "labelValues": label_values,
            "categories": attrs.get('categories', [])
        }

    async def _process_images_group(self, session: aiohttp.ClientSession, group_url: str, metadata: Dict[str, Any]):
        """Process the images group in SpatialData format."""
        logger.info(f"Processing images group at: {group_url}")
        
        # Look for common image names in SpatialData/Xenium format
        image_names = ['morphology_focus', 'morphology_mip', 'he', 'dapi']
        
        for image_name in image_names:
            try:
                image_url = urljoin(group_url, f'{image_name}/')
                zgroup_url = urljoin(image_url, '.zgroup')
                
                async with session.get(zgroup_url) as response:
                    if response.status == 200:
                        logger.info(f"Found image: {image_name}")
                        full_path = f"images/{image_name}"
                        metadata["datasetStructure"]["groups"].append(full_path)
                        
                        # Get image metadata
                        image_metadata = await self._extract_spatialdata_image_metadata(session, image_url, image_name)
                        metadata["images"][image_name] = image_metadata
                        
            except Exception as e:
                logger.debug(f"Image {image_name} not found: {e}")

    async def _process_points_group(self, session: aiohttp.ClientSession, group_url: str, metadata: Dict[str, Any]):
        """Process the points group in SpatialData format."""
        logger.info(f"Processing points group at: {group_url}")
        
        # Look for points data (typically transcripts in Xenium)
        points_names = ['transcripts']
        
        for points_name in points_names:
            try:
                points_url = urljoin(group_url, f'{points_name}/')
                zgroup_url = urljoin(points_url, '.zgroup')
                
                async with session.get(zgroup_url) as response:
                    if response.status == 200:
                        logger.info(f"Found points data: {points_name}")
                        full_path = f"points/{points_name}"
                        metadata["datasetStructure"]["groups"].append(full_path)
                        
                        # Get points metadata
                        points_metadata = await self._extract_spatialdata_points_metadata(session, points_url, points_name)
                        metadata["points"][points_name] = points_metadata
                        
            except Exception as e:
                logger.debug(f"Points data {points_name} not found: {e}")

    async def _process_shapes_group(self, session: aiohttp.ClientSession, group_url: str, metadata: Dict[str, Any]):
        """Process the shapes group in SpatialData format."""
        logger.info(f"Processing shapes group at: {group_url}")
        
        # Look for shape data (cell boundaries, circles, etc.)
        shapes_names = ['cell_boundaries', 'cell_circles', 'nuclei']
        
        for shape_name in shapes_names:
            try:
                shape_url = urljoin(group_url, f'{shape_name}/')
                zgroup_url = urljoin(shape_url, '.zgroup')
                
                async with session.get(zgroup_url) as response:
                    if response.status == 200:
                        logger.info(f"Found shape data: {shape_name}")
                        full_path = f"shapes/{shape_name}"
                        metadata["datasetStructure"]["groups"].append(full_path)
                        
                        # Get shapes metadata
                        shapes_metadata = await self._extract_spatialdata_shapes_metadata(session, shape_url, shape_name)
                        metadata["shapes"][shape_name] = shapes_metadata
                        
            except Exception as e:
                logger.debug(f"Shape data {shape_name} not found: {e}")

    async def _process_tables_group(self, session: aiohttp.ClientSession, group_url: str, metadata: Dict[str, Any]):
        """Process the tables group in SpatialData format."""
        logger.info(f"Processing tables group at: {group_url}")
        
        # Look for table data
        table_names = ['table']
        
        for table_name in table_names:
            try:
                table_url = urljoin(group_url, f'{table_name}/')
                zgroup_url = urljoin(table_url, '.zgroup')
                
                async with session.get(zgroup_url) as response:
                    if response.status == 200:
                        logger.info(f"Found table data: {table_name}")
                        full_path = f"tables/{table_name}"
                        metadata["datasetStructure"]["groups"].append(full_path)
                        
                        # Get table metadata - this is the most complex part
                        table_metadata = await self._extract_spatialdata_table_metadata(session, table_url, table_name)
                        metadata["tables"][table_name] = table_metadata
                        
            except Exception as e:
                logger.debug(f"Table data {table_name} not found: {e}")

    async def _extract_spatialdata_image_metadata(self, session: aiohttp.ClientSession, image_url: str, image_name: str) -> Dict[str, Any]:
        """Extract metadata from SpatialData image."""
        metadata = {
            "name": image_name,
            "type": "image",
            "scales": []
        }
        
        # Look for scale directories (0, 1, 2, etc.)
        for scale in range(5):  # Check scales 0-4
            try:
                scale_url = urljoin(image_url, f'{scale}/')
                zarray_url = urljoin(scale_url, '.zarray')
                
                async with session.get(zarray_url) as response:
                    if response.status == 200:
                        zarray_text = await response.text()
                        array_meta = json.loads(zarray_text)
                        
                        scale_info = {
                            "scale": scale,
                            "shape": array_meta.get('shape', []),
                            "dtype": array_meta.get('dtype', 'unknown'),
                            "chunks": array_meta.get('chunks', [])
                        }
                        metadata["scales"].append(scale_info)
                        
            except Exception as e:
                logger.debug(f"Scale {scale} not found for image {image_name}: {e}")
                break
        
        return metadata

    async def _extract_spatialdata_points_metadata(self, session: aiohttp.ClientSession, points_url: str, points_name: str) -> Dict[str, Any]:
        """Extract metadata from SpatialData points."""
        return {
            "name": points_name,
            "type": "points",
            "format": "parquet",
            "description": f"Point data for {points_name}"
        }

    async def _extract_spatialdata_shapes_metadata(self, session: aiohttp.ClientSession, shapes_url: str, shape_name: str) -> Dict[str, Any]:
        """Extract metadata from SpatialData shapes."""
        return {
            "name": shape_name,
            "type": "shapes",
            "format": "parquet",
            "description": f"Shape data for {shape_name}"
        }

    async def _extract_spatialdata_table_metadata(self, session: aiohttp.ClientSession, table_url: str, table_name: str) -> Dict[str, Any]:
        """Extract metadata from SpatialData table (most complex - AnnData-like structure)."""
        metadata = {
            "name": table_name,
            "type": "table",
            "components": {}
        }
        
        # Look for standard AnnData components
        components = ['X', 'obs', 'var', 'obsm', 'varm', 'layers', 'uns']
        
        for component in components:
            try:
                component_url = urljoin(table_url, f'{component}/')
                zgroup_url = urljoin(component_url, '.zgroup')
                
                async with session.get(zgroup_url) as response:
                    if response.status == 200:
                        logger.info(f"Found table component: {component}")
                        full_path = f"tables/{table_name}/{component}"
                        
                        if component == 'X':
                            # X matrix - sparse or dense
                            x_metadata = await self._extract_x_matrix_metadata(session, component_url)
                            metadata["components"][component] = x_metadata
                        elif component in ['obs', 'var']:
                            # Observation/variable metadata
                            metadata["components"][component] = {
                                "type": "dataframe",
                                "description": f"{component} annotations"
                            }
                        elif component in ['obsm', 'varm']:
                            # Multi-dimensional annotations
                            metadata["components"][component] = {
                                "type": "multidimensional", 
                                "description": f"{component} multidimensional data"
                            }
                        else:
                            metadata["components"][component] = {
                                "type": "group",
                                "description": f"{component} data"
                            }
                        
            except Exception as e:
                logger.debug(f"Table component {component} not found: {e}")
        
        return metadata

    async def _extract_x_matrix_metadata(self, session: aiohttp.ClientSession, x_url: str) -> Dict[str, Any]:
        """Extract metadata from X matrix (expression data)."""
        # Check for sparse matrix components
        sparse_components = ['data', 'indices', 'indptr']
        found_sparse = 0
        
        for component in sparse_components:
            try:
                comp_url = urljoin(x_url, f'{component}/')
                zarray_url = urljoin(comp_url, '.zarray')
                
                async with session.get(zarray_url) as response:
                    if response.status == 200:
                        found_sparse += 1
                        
            except Exception:
                pass
        
        if found_sparse == 3:
            return {
                "type": "sparse_matrix",
                "format": "csr", 
                "description": "Sparse expression matrix"
            }
        else:
            # Try dense matrix
            try:
                zarray_url = urljoin(x_url, '.zarray')
                async with session.get(zarray_url) as response:
                    if response.status == 200:
                        return {
                            "type": "dense_matrix",
                            "description": "Dense expression matrix"
                        }
            except Exception:
                pass
        
        return {
            "type": "unknown",
            "description": "Expression matrix"
        }


# Global instance for reuse
_extractor_instance = None

def get_zarr_extractor() -> ZarrMetadataExtractor:
    """Get or create global ZarrMetadataExtractor instance"""
    global _extractor_instance
    if _extractor_instance is None:
        _extractor_instance = ZarrMetadataExtractor()
    return _extractor_instance