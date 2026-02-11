import argparse
import os
from dataclasses import dataclass
from typing import TYPE_CHECKING, Optional
import numpy as np
from ..serverlite import serve_project
# nb, main spatialdata import should happen lazily
# so user doesn't have to wait and see lots of scary unrelated output if they input bad arguments.
if TYPE_CHECKING:
    # import spatialdata as sd
    from anndata import AnnData
    from spatialdata.transformations import BaseTransformation
    from spatialdata import SpatialData
    from spatialdata.models import SpatialElement
    from mdvtools.mdvproject import MDVProject
    from geopandas import GeoDataFrame

@dataclass
class SpatialDataConversionArgs:
    spatialdata_path: str
    output_folder: str
    temp_folder: str
    preserve_existing: bool = False
    output_geojson: bool = False
    serve: bool = False
    link: bool = False
    point_transform: str = "auto"

def _process_sdata_path(sdata_path: str, conversion_args: "SpatialDataConversionArgs"):
    """Processes a single SpatialData object path."""
    # imports need to be here for the separate process
    from mdvtools.spatial.conversion import _try_read_zarr, _resolve_regions_for_table
    from mdvtools.spatial.mermaid import sdata_to_mermaid
    import os

    sdata_name = os.path.basename(sdata_path)
    sdata = _try_read_zarr(sdata_path)
    if sdata is None:
        return None
    print(f"## SpatialData object representation:\n\n```\n{sdata}\n```\n")
    print(f"## Mermaid diagram:\n\n{sdata_to_mermaid(sdata)}\n")
    adata_objects = []
    all_regions = {}
    for table_name, adata in sdata.tables.items():
        _resolve_regions_for_table(sdata, table_name, sdata_name, conversion_args)
        adata.obs["spatialdata_path"] = sdata_name
        adata.obs["table_name"] = table_name
        adata_objects.append(adata)
        if "regions" in adata.uns.get("mdv", {}):
            all_regions.update(adata.uns["mdv"]["regions"])

    return sdata_name, sdata, adata_objects, all_regions


REGION_FIELD = "spatial_region"
@dataclass
class ImageEntry:
    """
    A single image that is associated with a region in an AnnData table.
    Local helper, not part of the public API or set in stone.
    """
    path: str
    transform_to_image: "BaseTransformation"
    extent_px: tuple[int, int]
    is_primary: bool
    region_id: str

def add_readme_to_project(mdv: "MDVProject", adata: Optional["AnnData"], conversion_args: Optional["SpatialDataConversionArgs"] = None, initial_markdown = ""):
    """
    Add a README.md file to the project.

    We may consider adding more configuration options etc.
    """
    from mdvtools.spatial.mermaid import sdata_to_mermaid
    from mdvtools.markdown_utils import create_project_markdown
    from mdvtools.charts.text_box_plot import TextBox
    from mdvtools.build_info import get_build_info_markdown
    markdown = initial_markdown
    markdown += "# SpatialData conversion information:\n\n"
    
    # Add build information
    markdown += get_build_info_markdown()
    
    # Add conversion arguments summary
    if conversion_args is not None:
        markdown += "### Conversion arguments:\n\n"
        markdown += f"- **Point transform mode**: `{conversion_args.point_transform}`\n"
        markdown += f"- **Preserve existing**: {conversion_args.preserve_existing}\n"
        markdown += f"- **Output GeoJSON**: {conversion_args.output_geojson}\n"
        markdown += f"- **Link spatial data**: {conversion_args.link}\n"
        markdown += f"- **SpatialData path**: `{conversion_args.spatialdata_path}`\n"
        markdown += "\n"
    
    markdown += "## SpatialData objects:\n\n"
    spatial_dir = os.path.join(mdv.dir, "spatial")
    
    for sdata_name in os.listdir(spatial_dir):
        sdata_path = os.path.join(spatial_dir, sdata_name)
        if not os.path.exists(sdata_path):
            continue
        sdata = _try_read_zarr(sdata_path)
        if sdata is None:
            continue
        sdata_md = f"## {sdata_name}:\n\n```\n{sdata}\n```"
        sdata_md += f"\n\n### Mermaid diagram:\n\n{sdata_to_mermaid(sdata)}\n\n"
        markdown += sdata_md
    
    # Add point transform strategy section
    if adata is not None and "mdv" in adata.uns and "point_transform" in adata.uns["mdv"]:
        markdown += "\n\n---\n\n## Point transform strategy:\n\n"
        markdown += "The following point transform strategies were used for coordinate conversion:\n\n"
        
        point_transform_info = adata.uns["mdv"]["point_transform"]
        if isinstance(point_transform_info, dict):
            for region_id, transform_meta in point_transform_info.items():
                xenium_note = " (Xenium detected)" if transform_meta.get("xenium_detected", False) else ""
                markdown += f"- **{region_id}**:\n"
                markdown += f"  - Mode: `{transform_meta['mode']}`{xenium_note}\n"
                markdown += f"  - Source element: `{transform_meta['source_element']}`\n"
                markdown += f"  - Target element: `{transform_meta['target_element']}`\n"
                markdown += f"  - Transform type: `{transform_meta['transform_type']}`\n"
        
        markdown += "\n"
    
    markdown += "\n\n---\n\n## Project DataSource summary:\n\n"
    markdown += create_project_markdown(mdv, False)
    # todo - add anndata markdown here
    # if adata is not None:
    #     markdown += f"\n{create_anndata_markdown(adata)}"
    with open(os.path.join(mdv.dir, "README.md"), "w") as f:
        f.write(markdown)
    print(f"README.md written to {os.path.join(mdv.dir, 'README.md')}")
    textbox = TextBox(title="SpatialData conversion information", param=[], size=[792, 472], position=[10, 10])
    textbox.set_text(markdown)
    mdv.set_view("Data summary", {"initialCharts": {"cells": [textbox.plot_data]}})

def _shape_to_geojson(shape: "GeoDataFrame", transform_to_image: "BaseTransformation") -> str:
    """
    Convert a spatial element to a GeoJSON object, with the coordinates transformed by the provided transformation.
    """
    axes = ("x", "y")
    T = transform_to_image.to_affine_matrix(input_axes=axes, output_axes=axes)
    # The affine matrix T is a 3x3 matrix.
    # The affine_transform function from geopandas takes a 6-element tuple or list of coefficients:
    # [a, b, d, e, xoff, yoff]
    # which correspond to the transformation matrix:
    # | a  b  xoff |
    # | d  e  yoff |
    # | 0  0  1    |
    # We extract these coefficients from the matrix T.
    matrix = [T[0, 0], T[0, 1], T[1, 0], T[1, 1], T[0, 2], T[1, 2]]
    transformed_shape = shape.affine_transform(matrix)
    return f"{transformed_shape.to_json()}"

def _get_transformation_between_coordinate_systems_safe(
    sdata: "SpatialData",
    source_coordinate_system: "SpatialElement | str",
    target_coordinate_system: "SpatialElement | str",
) -> Optional["BaseTransformation"]:
    """
    Safely get transformation between coordinate systems, handling multiple equal paths.
    
    When multiple equal-length shortest paths are found, this function automatically
    selects the first path and logs a warning, rather than raising an error.
    
    Args:
        sdata: SpatialData object
        source_coordinate_system: Source coordinate system (element or string)
        target_coordinate_system: Target coordinate system (element or string)
        
    Returns:
        Transformation between the coordinate systems, or None if no path exists
    """
    from spatialdata.transformations import get_transformation_between_coordinate_systems
    from spatialdata.transformations import get_transformation, Identity, Sequence
    from spatialdata.models._utils import has_type_spatial_element
    import networkx as nx
    import contextlib
    
    try:
        return get_transformation_between_coordinate_systems(
            sdata, source_coordinate_system, target_coordinate_system
        )
    except RuntimeError as e:
        error_msg = str(e)
        if "Multiple equal paths found" not in error_msg:
            # Re-raise if it's a different error
            raise
        
        # Log warning about auto-selecting path
        print(
            f"WARNING: Multiple equal transformation paths found. "
            f"Automatically selecting the first path. Error details: {error_msg}"
        )
        
        # Manually build graph and select first shortest path
        # Replicate _build_transformations_graph logic
        g = nx.DiGraph()
        gen = sdata._gen_spatial_element_values()
        for cs in sdata.coordinate_systems:
            g.add_node(cs)
        for e in gen:
            g.add_node(id(e))
            transformations = get_transformation(e, get_all=True)
            if not isinstance(transformations, dict):
                raise ValueError("Expected transformations to be a dict")
            for cs, t in transformations.items():
                g.add_edge(id(e), cs, transformation=t)
                with contextlib.suppress(np.linalg.LinAlgError):
                    g.add_edge(cs, id(e), transformation=t.inverse())
        
        # Determine source and target nodes
        if has_type_spatial_element(source_coordinate_system):
            src_node = id(source_coordinate_system)
        else:
            assert isinstance(source_coordinate_system, str)
            src_node = source_coordinate_system
        
        if has_type_spatial_element(target_coordinate_system):
            tgt_node = id(target_coordinate_system)
        else:
            assert isinstance(target_coordinate_system, str)
            tgt_node = target_coordinate_system
        
        # Check for identity case
        if src_node == tgt_node:
            return Identity()
        
        try:
            path = nx.shortest_path(g, source=src_node, target=tgt_node)
            
            # Build transformation sequence from path
            transformations_list = [
                g[path[i]][path[i + 1]]["transformation"] 
                for i in range(len(path) - 1)
            ]
            sequence = Sequence(transformations_list)
            return sequence
        
        # - not expected as we call this with nodes previously established to be connected
        except nx.NetworkXNoPath:
            # No path exists between source and target 
            return None
        except Exception as e:
            # Re-raise all other exceptions as unexpected errors
            raise RuntimeError(
                f"Unexpected error while computing transformation path: {e}"
            ) from e


def _transform_table_coordinates(adata: "AnnData", region_to_image: dict[str, ImageEntry]):
    """
    Transform the coordinates of an AnnData table to the coordinates of the images associated with the regions.
    """
    from spatialdata.models import get_table_keys
    from spatialdata.transformations import Identity
    if "spatial" not in adata.obsm:
        # todo: proper support for non-spatial tables, add tests.
        print("WARNING: No spatial coordinates found in obsm['spatial']")
        print("We should be able to handle this case, but it is not supported or tested yet, results may be undesired.")
        return adata
    # todo: support non-2d cases...
    axes = ("x", "y")
    # Convert coordinates to homogeneous coordinates (add 1s for translation)
    coords_homogeneous = np.column_stack([adata.obsm["spatial"], np.ones(adata.obsm["spatial"].shape[0], dtype=float)])
    _r, region_element, _instance_key = get_table_keys(adata)

    transformed_coords = np.full_like(adata.obsm["spatial"], fill_value=np.nan)
    region_ids = np.empty(adata.n_obs, dtype=object)

    for r, entry in region_to_image.items():
        if region_element is None:
            mask = np.ones(adata.n_obs, dtype=bool)
        else:
            mask = (adata.obs[region_element] == r).values
        if not np.any(mask):
            continue

        # Handle Identity transform specially (just copy coordinates)
        if isinstance(entry.transform_to_image, Identity):
            transformed_coords[mask, :] = adata.obsm["spatial"][mask, :]
        else:
            T = entry.transform_to_image.to_affine_matrix(input_axes=axes, output_axes=axes)
            transformed_coords_homogeneous = coords_homogeneous[mask] @ T.T
            transformed_coords[mask, :] = transformed_coords_homogeneous[:, :-1]
        
        region_ids[mask] = entry.region_id

    adata.obs["x"] = transformed_coords[:, 0]
    adata.obs["y"] = transformed_coords[:, 1]
    adata.obs[REGION_FIELD] = region_ids
    return adata

def _get_transform_keys(e: "SpatialElement") -> list[str]:
    """
    Returns a list of coordinate system names that are defined for the given element.
    """
    from spatialdata.transformations import get_transformation
    transformations = get_transformation(e, get_all=True)
    if not isinstance(transformations, dict):
        raise ValueError("This should be unreachable, get_transformation with get_all=True should always return a dict")
    return list(transformations.keys())


def _is_xenium_like(sdata: "SpatialData") -> bool:
    """
    Check if SpatialData object appears to be Xenium data.
    
    Detection is based on presence of 'morphology_focus' image, which is
    a robust indicator of Xenium data.
    """
    return "morphology_focus" in sdata.images


def _find_xenium_shape(sdata: "SpatialData") -> Optional[tuple[str, "SpatialElement"]]:
    """
    Find a Xenium shape element (cell_circles, cell_boundaries, or nucleus_boundaries).
    
    Returns:
        Tuple of (name, element) if found, None otherwise.
    """
    xenium_shapes = ["cell_circles", "cell_boundaries", "nucleus_boundaries"]
    for shape_name in xenium_shapes:
        if shape_name in sdata.shapes:
            return (shape_name, sdata.shapes[shape_name])
    return None


def _describe_transform_type(transform: "BaseTransformation") -> str:
    """
    Get a human-readable description of a transformation.
    
    Args:
        transform: The transformation to describe
        
    Returns:
        String description of the transform type
    """
    from spatialdata.transformations import Identity, Scale
    
    if isinstance(transform, Identity):
        return "Identity"
    elif isinstance(transform, Scale):
        scale = transform.scale
        if len(scale) == 2:
            if scale[0] == scale[1]:
                return f"Scale({scale[0]:.2f})"
            else:
                return f"Scale({scale[0]:.2f}, {scale[1]:.2f})"
    return type(transform).__name__


def _get_xenium_transform(
    sdata: "SpatialData",
    img_obj: "SpatialElement",
    img_name: str,
    sdata_name: str,
    require: bool = False
) -> Optional[tuple["BaseTransformation", str]]:
    """
    Get transform from Xenium shape to image, if Xenium data is detected.
    
    Args:
        sdata: SpatialData object
        img_obj: Image element to transform to
        img_name: Name of the image element
        sdata_name: Optional name for error messages
        require: If True, raise errors if Xenium not found; if False, return None
        
    Returns:
        Tuple of (transform, shape_name) if found, None otherwise
    """
    if not _is_xenium_like(sdata):
        if require:
            context = f" in '{sdata_name}'" if sdata_name else ""
            raise ValueError(
                f"Xenium mode specified but 'morphology_focus' image not found{context}. "
                f"This mode requires Xenium data structure."
            )
        return None
    
    xenium_shape = _find_xenium_shape(sdata)
    if xenium_shape is None:
        if require:
            context = f" in '{sdata_name}'" if sdata_name else ""
            raise ValueError(
                f"Xenium mode specified but no Xenium shapes (cell_circles, cell_boundaries, "
                f"nucleus_boundaries) found{context}."
            )
        return None
    
    shape_name, shape_element = xenium_shape
    transform = _get_transformation_between_coordinate_systems_safe(sdata, shape_element, img_obj)
    if transform is None:
        if require:
            context = f" in '{sdata_name}'" if sdata_name else ""
            raise ValueError(
                f"Could not compute transform from Xenium shape '{shape_name}' to image '{img_name}'{context}."
            )
        return None
    
    return (transform, shape_name)


def _choose_point_transform(
    sdata: "SpatialData",
    annotated_element: "SpatialElement",
    annotated_name: str,
    img_obj: "SpatialElement",
    img_name: str,
    mode: str,
    sdata_name: str = ""
) -> tuple[Optional["BaseTransformation"], dict]:
    """
    Choose the appropriate point transform based on the specified mode.
    
    Args:
        sdata: SpatialData object
        annotated_element: The element annotated by the table
        annotated_name: Name of the annotated element
        img_obj: Image element to transform to
        img_name: Name of the image element
        mode: Transform mode ('image', 'auto', 'xenium', 'identity', 'annotated-element')
        sdata_name: Optional name of the SpatialData object (for error messages)
        
    Returns:
        Tuple of (transform, metadata_dict) where metadata contains:
        - mode: the strategy used
        - source_element: name of element used
        - target_element: name of image used
        - transform_type: brief description
        - xenium_detected: boolean if auto-mode detected Xenium
    """
    metadata = {
        "mode": mode,
        "source_element": annotated_name,
        "target_element": img_name,
        "transform_type": None,
        "xenium_detected": False
    }
    
    # Handle identity mode
    if mode == "identity":
        metadata["transform_type"] = "Identity"
        return None, metadata
    
    # Handle xenium mode (explicit)
    if mode == "xenium":
        result = _get_xenium_transform(sdata, img_obj, img_name, sdata_name, require=True)
        if result is None:
            raise ValueError("Unexpected: _get_xenium_transform returned None with require=True")
        transform, shape_name = result
        metadata.update({
            "source_element": shape_name,
            "transform_type": _describe_transform_type(transform),
            "xenium_detected": True
        })
        return transform, metadata
    
    # Handle auto mode (detect Xenium, fallback to image)
    if mode == "auto":
        result = _get_xenium_transform(sdata, img_obj, img_name, sdata_name, require=False)
        if result is not None:
            transform, shape_name = result
            metadata.update({
                "mode": "auto (xenium)",
                "source_element": shape_name,
                "transform_type": _describe_transform_type(transform),
                "xenium_detected": True
            })
            return transform, metadata
        # Fallback to image mode
        mode = "image"
        metadata["mode"] = "auto (image)"
    
    # Handle annotated-element mode
    if mode == "annotated-element":
        # Note: This transform is from the element's coordinate system, not to image
        # We still need to get the transform to the image
        transform = _get_transformation_between_coordinate_systems_safe(sdata, annotated_element, img_obj)
        if transform is None:
            metadata["transform_type"] = "None (no transform to image found)"
        else:
            metadata["transform_type"] = _describe_transform_type(transform)
        return transform, metadata
    
    # Handle image mode (default)
    if mode == "image":
        transform = _get_transformation_between_coordinate_systems_safe(sdata, annotated_element, img_obj)
        if transform is None:
            metadata["transform_type"] = "None (no transform found)"
        else:
            metadata["transform_type"] = _describe_transform_type(transform)
        return transform, metadata
    
    raise ValueError(f"Unknown point_transform mode: '{mode}'")

# ---- region resolution ----
def _resolve_regions_for_table(sdata: "SpatialData", table_name: str, sdata_name: str, conversion_args: "SpatialDataConversionArgs"):
    """
    Internal processing of a table with various side-effects on the provided adata object:
    - the annotated elements for each region (which may be shapes etc) is used to establish an 'image' that should be appropriate for the region.
    - a new set of spatial coordinates are added to obs[x, y], with the obs['spatial'] coordinates transformed to the image coordinates.
      (this is so that current MDV - which isn't able to process transformation from spatialdata - can display the data. The original data is preserved)
    - metadata is added to uns["mdv"]["regions"] in a form that resembles what will be needed for mdv regions
      this is used later to set the region metadata for the mdv project.
    If conversion_args.output_geojson is True, geojson should be saved to the output folder as <mdv.dir>/spatial/region_id.geojson
    """
    from spatialdata.transformations import get_transformation
    from spatialdata.models import get_table_keys
    from anndata import AnnData
    adata = sdata.tables[table_name]
    if not isinstance(adata, AnnData):
        # if you hit this, it is likely because an invalid table_name was passed to this function.
        raise ValueError(f"Invalid table_name? '{table_name}' is not a table in '{sdata_name}':\n{sdata}")

    adata.uns.setdefault("mdv", {})
    # Early return for non-spatial tables
    if "spatial" not in adata.obsm:
        print(f"INFO: Skipping region resolution for non-spatial table '{table_name}' in '{sdata_name}'")
        # Mark as non-spatial for downstream handling
        adata.uns["mdv"]["is_spatial"] = False
        return adata

    # Mark as spatial
    adata.uns["mdv"]["is_spatial"] = True

    
    region, _obs_region_col, _instance_key = get_table_keys(adata)

    regions = region if isinstance(region, list) else [region]
    if not regions:
        raise ValueError(f"No region for table '{table_name}'")

    all_regions: dict[str, dict] = {}

    region_to_image: dict[str, ImageEntry] = {}
    transform_provenance: dict[str, dict] = {}  # region_id -> provenance metadata
    
    for r in regions:
        annotated = sdata[r]  # shapes/points/labels/image ... probably shapes, but we need to find an appropriate image.
        if isinstance(annotated, AnnData):
            raise ValueError(f"This should be unreachable - '{r}' is not a table in '{sdata_name}':\n{sdata}")
        transformations = get_transformation(annotated, get_all=True)
        if not isinstance(transformations, dict):
            raise ValueError("This should be unreachable, get_transformation with get_all=True should always return a dict")
        cs_names = list(transformations.keys())
        if not cs_names:
            raise ValueError(f"Annotated element '{r}' has no coordinate system definitions")

        # Heuristic: prefer a CS literally named like the region, else the first
        cs = r if r in cs_names else cs_names[0]

        # images that are likely compatible
        img_candidates = [img for img in sdata.images.items() if cs in _get_transform_keys(img[1])]

        # Sort to prioritize morphology_focus images (especially for xenium datasets)
        # note - there is some dead code below that would interact badly with this if it was actually doing anything
        img_candidates.sort(key=lambda img: (0 if 'morphology_focus' in str(img[0]) else 1, img[0]))

        img_entries: list[ImageEntry] = []
        best_idx = 0
        best_transform_metadata = None
        
        for i, (img_path, img_obj) in enumerate(img_candidates):
            # Use the chosen transform strategy
            T, transform_metadata = _choose_point_transform(
                sdata, annotated, r, img_obj, img_path,
                conversion_args.point_transform, sdata_name
            )
            if T is None and conversion_args.point_transform != "identity":
                continue
            
            ## TODO fix or simplify non-functional size-heuristic...
            # currently not expecting to do much with wh, and anticipating that this will be less relevant
            # once we use spatialdata.js
            
            # try to extract pixel size / extent if available
            # nb, this is wrong, it was written by an LLM.
            # not really used anyway so just putting some arbitrary defaults.
            # extent_px = getattr(img_obj, "shape", None)  # (C?, Y, X) or (Y, X)
            # if extent_px is not None:
            #     wh = (extent_px[-1], extent_px[-2])
            # else:
            #     print(f"WARNING: No extent found for image '{img_path}' - using default 1000x1000")
            #     print("This is only used for legacy reasons anyway and is unlikely to impact functionality.")
            #     wh = (1000, 1000)
            wh = (1000, 1000)

            # For identity mode, we still need a transform object for ImageEntry
            # but we'll handle it specially in _transform_table_coordinates
            if T is None:
                from spatialdata.transformations import Identity
                T = Identity()

            img_entries.append(ImageEntry(
                path=img_path,
                transform_to_image=T,
                extent_px=wh,
                is_primary=False,
                region_id=f"{sdata_name}_{r}"
            ))
            
            # Store metadata for the first viable transform (we'll use best_idx later)
            if best_transform_metadata is None:
                best_transform_metadata = transform_metadata
            
            # primary selection heuristic: first viable with largest area
            prev = img_entries[best_idx].extent_px
            if wh and prev and (wh[0]*wh[1] > prev[0]*prev[1]):
                best_idx = len(img_entries)-1
                best_transform_metadata = transform_metadata


        if not img_entries:
            # raise ValueError(f"No images found for region '{r}' in '{sdata_name}'")
            print(f"WARNING: No images found for region '{r}' in '{sdata_name}'")
            continue
        
        for j, entry in enumerate(img_entries):
            entry.is_primary = (j == best_idx)

        # todo - also add images that aren't the primary image...
        best_img = img_entries[best_idx]
        region_id = best_img.region_id
        
        # Store provenance metadata
        if best_transform_metadata:
            transform_provenance[region_id] = best_transform_metadata
            # Log the transform decision
            xenium_note = " (Xenium detected)" if best_transform_metadata.get("xenium_detected", False) else ""
            print(
                f"INFO: Table '{table_name}' region '{r}' using point_transform mode '{best_transform_metadata['mode']}'"
                f"{xenium_note}: {best_transform_metadata['source_element']} -> {best_transform_metadata['target_element']} "
                f"({best_transform_metadata['transform_type']})"
            )
        all_regions[region_id] = {
            "spatial": {
                "coordinate_system": cs,
                "file": sdata_name
            },
            "viv_image": {
                # relative to the "avivator" base set later
                "file": f"{sdata_name}/images/{best_img.path}"
            },
            "roi": {
                "min_x": 0,
                "min_y": 0,
                "max_x": best_img.extent_px[0],
                "max_y": best_img.extent_px[1],
            },
            "images": {},
        }

        # do we want to save geojson for this region?
        # nb this will be deprecated once we have spatialdata.js layers with shapes.
        if conversion_args.output_geojson:
            from geopandas import GeoDataFrame
            # xenium hack... we want cell_boundaries, not cell_circles which is the annotated element...
            if "cell_boundaries" in sdata.shapes:
                print(f"Xenium hack... using cell_boundaries for geojson output for region '{r}' in '{sdata_name}'")
                annotated = sdata["cell_boundaries"]
            
            if isinstance(annotated, GeoDataFrame):
                geojson = _shape_to_geojson(annotated, best_img.transform_to_image)
                region_id = best_img.region_id
                # the name may be misleading in xenium case, but avoids other potential naming conflicts.
                name = f"{region_id}.geo.json"
                path = os.path.join(conversion_args.temp_folder, name)
                all_regions[region_id]["json"] = f"images/{name}"
                with open(path, "w") as f:
                    f.write(geojson)
                    print(f"Wrote geojson for region '{region_id}' to {path}")
            else:
                print(f"WARNING: No geojson output for region '{r}' in '{sdata_name}' because it is not a GeoDataFrame")
        
        region_to_image[r] = best_img

    # transform coordinates to image coordinates for each row in the table
    _transform_table_coordinates(adata, region_to_image)
    
    # store in uns for this table
    adata.uns.setdefault("mdv", {})
    adata.uns["mdv"].setdefault("regions", all_regions)
    # Store transform provenance
    adata.uns["mdv"]["point_transform"] = transform_provenance




def _try_read_zarr(path: str):# -> sd.SpatialData | None:
    """
    Attempts to read arbitrary path string as a SpatialData object, falling back to None if it fails.
    """
    try:
        import spatialdata as sd
        return sd.read_zarr(path)
    except Exception as e:
        print(f"Warning: Failed to read SpatialData object from {path}: '{e}'")
        return None

def _set_default_image_view(mdv: "MDVProject"):
    """
    Uses a template view for the default view of the spatial data.
    In future it would be good to make this user-configurable.
    """
    import json
    import os

    # find a region with a viv_image and use key from that...
    for region_name, region_data in mdv.get_datasource_metadata("cells")["regions"]["all_regions"].items():
        if "viv_image" in region_data:
            break
    else:
        raise ValueError("No region with a viv_image found")

    print(f"Using region '{region_name}' for default view")
    with open(os.path.join(os.path.dirname(__file__), "spatial_view_template.json"), "r") as f:
        view_str = f.read().replace("<SPATIAL_REGION_NAME>", region_name)
        mdv.set_view("default", json.loads(view_str), True)

def convert_spatialdata_to_mdv(args: SpatialDataConversionArgs):
    """
    Convert all SpatialData objects in a folder to MDV format.
    Args:
        spatialdata_path: Path to the folder containing the SpatialData objects
        output_folder: Path to the output folder for the MDV project
        preserve_existing: Whether to preserve existing project data in output folder
        serve: Whether to serve the project after conversion
    Returns:
        MDVProject: The MDV project object
    Raises:
        ValueError: If the SpatialData objects cannot be found or the project cannot be created
        Exception: For other unexpected errors during conversion
    """
    # imports can be slow, so doing them here rather than at the top of the file
    from anndata import concat as ad_concat
    from concurrent.futures import ProcessPoolExecutor

    # from mdvtools.spatial.spatial_conversion import convert_spatialdata_to_mdv
    from mdvtools.conversions import convert_scanpy_to_mdv
    from mdvtools.markdown_utils import create_project_markdown
    from mdvtools.build_info import get_build_info
    
    # Print version banner if build info is available
    build_info = get_build_info()
    if build_info["source"] != "unknown" and build_info["git_commit_hash"]:
        commit_short = build_info["git_commit_hash"][:8]
        date_str = build_info["git_commit_date"] or build_info["build_date"] or "unknown"
        print(f"MDV conversion (commit: {commit_short}, date: {date_str})")
    

    # we could do a nicer glob thing here, but I don't want to test that right now.
    sdata_paths = [
        os.path.join(args.spatialdata_path, f)
        for f in os.listdir(args.spatialdata_path)
    ]
    # sdata_paths = [f for f in sdata_paths if f.endswith(".zarr")]
    sdata_paths = sorted(sdata_paths)
    if len(sdata_paths) == 0:
        # we haven't actually yet determined whether any of the paths are really sdata...
        # this happens when the folder is totally empty.
        raise ValueError(f"No SpatialData objects found in the folder '{args.spatialdata_path}'")
    sdata_objects: dict[str, SpatialData] = {}
    adata_objects: list[AnnData] = []
    all_regions: dict[str, dict] = {}
    names: set[str] = set()
    with ProcessPoolExecutor() as executor:
        results = executor.map(_process_sdata_path, sdata_paths, [args] * len(sdata_paths))

    for result in results:
        if result is None:
            continue
        
        sdata_name, sdata, adatas, regions = result
        if sdata_name in names:
            raise ValueError(
                f"SpatialData object with name '{sdata_name}' already processed. Please ensure names are unique."
            )
        
        names.add(sdata_name)
        sdata_objects[sdata_name] = sdata
        adata_objects.extend(adatas)
        all_regions.update(regions)

    if len(adata_objects) == 0:
        raise ValueError("No tables found in any SpatialData objects - this is not yet supported.")
    
    # Check if we have at least one spatial table - maybe this is unnecessary noise?
    spatial_tables = [ad for ad in adata_objects if ad.uns.get("mdv", {}).get("is_spatial", False)]

    if len(spatial_tables) == 0:
        print("WARNING: No spatial tables found. Skipping spatial-specific metadata and view setup.")
        # Still merge and convert, but skip spatial operations
        merged_adata = ad_concat(adata_objects, index_unique="_")
        mdv = convert_scanpy_to_mdv(
            args.output_folder, merged_adata, delete_existing=not args.preserve_existing
        )
        # ... write sdata objects but skip region metadata and _set_default_image_view
        return mdv

    # Proceed with normal spatial workflow if we have spatial tables
    if not all_regions:
        raise ValueError("Spatial tables found but no regions could be resolved - this indicates a problem with the data")

    ## todo - try to make sure we have sparse CSC matrices if possible.
    merged_adata = ad_concat(adata_objects, index_unique="_")
    mdv = convert_scanpy_to_mdv(
        args.output_folder, merged_adata, delete_existing=not args.preserve_existing
    )
    if args.link:
        print(f"Linking spatialdata object path '{args.spatialdata_path}' to '{mdv.dir}/spatial'")
        os.symlink(
            os.path.abspath(args.spatialdata_path), 
            os.path.join(mdv.dir, "spatial"), target_is_directory=True
        )
    else:
        print(f"Copying spatialdata object path '{args.spatialdata_path}' to '{mdv.dir}/spatial'")
        os.makedirs(
            f"{mdv.dir}/spatial", exist_ok=True
        )  # pretty sure sdata.write will do this anyway
        for sdata_path, sdata in sdata_objects.items():
            sdata_name = os.path.basename(sdata_path)
            # pay attention to format specifier here?
            sdata.write(os.path.join(mdv.dir, "spatial", sdata_name))
    # move contents of temp_folder to output_folder/spatial
    # nb would probably rather avoid non spatialdata in spatial folder, 
    # especially if linking (should avoid side-effects in linked folder)
    if args.output_geojson:
        import shutil
        dest_dir = os.path.join(mdv.dir, "images")
        os.makedirs(dest_dir, exist_ok=True)
        for f in os.listdir(args.temp_folder):
            shutil.move(os.path.join(args.temp_folder, f), os.path.join(dest_dir, f))

    mdv.set_editable()

    # these methods won't do what we actually want... this should be addressed.
    # mdv.set_region_data("cells", all_regions,
    #                     region_field="coordinate_system",
    #                     default_color="x", # todo: figure out a better way to determine this.
    #                     position_fields=["x", "y"])
    #
    # mdv.add_viv_viewer("cells", [])
    # viv_data = {
    #     "path": f"spatial/{sdata_path}/images/{main_image}",
    # }
    # mdv.add_viv_images("cells", viv_data, link_images=False)

    # instead, we will set the region metadata directly.
    cells_md = mdv.get_datasource_metadata("cells")
    cells_md["regions"] = {
        # non-2d cases are one thing... but mixtures are another.
        # this metadata format won't work for that, we expect to move more towards spatialdata.js.
        "position_fields": ["x", "y"],
        "region_field": REGION_FIELD,
        "default_color": "x",  # todo: figure out a better way to determine this.
        "scale_unit": "µm",
        # scale will be baked during conversion from spatialdata to mdv.
        # in future we will use spatialdata.js to handle transforms at runtime.
        "scale": 1.0,
        "avivator": {
            "default_channels": [],
            "base_url": "spatial/",
        },
        "all_regions": all_regions,
    }
    mdv.set_datasource_metadata(cells_md)
    _set_default_image_view(mdv)
    print(f"## Merged AnnData object representation:\n\n```\n{merged_adata}\n```\n")
    print(f"## Project markdown:\n\n{create_project_markdown(mdv, False)}\n\n---")
    add_readme_to_project(mdv, merged_adata, args)
    if args.serve:
        print(f"Serving project at {args.output_folder}")
        serve_project(mdv)
    else:
        print(f"Project saved to {args.output_folder}")

    return mdv
if __name__ == "__main__":
    # take an array of spatialdata objects from a folder and convert them to mdv.
    parser = argparse.ArgumentParser(description="Convert SpatialData to MDV format")
    parser.add_argument("spatialdata_path", type=str, help="Path to SpatialData data")
    parser.add_argument("output_folder", type=str, help="Output folder for MDV project")
    parser.add_argument("--link", action="store_true", help="Symlink to the original SpatialData objects")
    parser.add_argument("--preserve-existing", action="store_true", help="Preserve existing project data")
    parser.add_argument("--output_geojson", action="store_true", help="Output geojson for each region (this feature to be deprecated in favour of spatialdata.js layers with shapes)")
    parser.add_argument("--serve", action="store_true", help="Serve the project after conversion")
    parser.add_argument(
        "--point-transform",
        type=str,
        default="auto",
        choices=["image", "auto", "xenium", "identity", "annotated-element"],
        help=(
            "Strategy for transforming point coordinates in tables to image coordinates. "
            "Options: 'image' (use annotated element to image transform), "
            "'auto' (detect Xenium and use shape-based transform, else fallback to 'image'), "
            "'xenium' (explicitly treat as Xenium, error if not found), "
            "'identity' (no transform, copy coordinates directly), "
            "'annotated-element' (use transformation from annotated element itself)."
        )
    )
    args = parser.parse_args()
    print(f"Converting SpatialData from '{args.spatialdata_path}' to {args.output_folder}")
    print(f"Preserving existing project data: {args.preserve_existing}")
    print(f"Serving the project after conversion: {args.serve}")
    print(f"Converting {len(os.listdir(args.spatialdata_path))} potential SpatialData objects")
    import tempfile
    with tempfile.TemporaryDirectory() as temp_folder:
        args.temp_folder = temp_folder
        convert_spatialdata_to_mdv(args) # type: ignore argparse->SpatialDataConversionArgs
        
