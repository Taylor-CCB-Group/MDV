import argparse
import os
from dataclasses import dataclass
from typing import TYPE_CHECKING
import numpy as np
# nb, main spatialdata import should happen lazily
# so user doesn't have to wait and see lots of scary unrelated output if they input bad arguments.
if TYPE_CHECKING:
    # import spatialdata as sd
    from anndata import AnnData
    from spatialdata.transformations import BaseTransformation
    from spatialdata import SpatialData
    from spatialdata.models import SpatialElement

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



def _transform_table_coordinates(adata: "AnnData", region_to_image: dict[str, ImageEntry]):
    """
    Transform the coordinates of an AnnData table to the coordinates of the images associated with the regions.
    """
    from spatialdata.models import get_table_keys
    if "spatial" not in adata.obsm:
        # todo: proper support for non-spatial tables, add tests.
        print(f"WARNING: No spatial coordinates found in obsm['spatial']")
        print(f"We should be able to handle this case, but it is not supported or tested yet, results may be undesired.")
        return adata
    # todo: support non-2d cases...
    axes = ("x", "y")
    # Convert coordinates to homogeneous coordinates (add 1s for translation)
    coords_homogeneous = np.column_stack([adata.obsm["spatial"], np.ones(adata.obsm["spatial"].shape[0])])
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

        T = entry.transform_to_image.to_affine_matrix(input_axes=axes, output_axes=axes)
        transformed_coords_homogeneous = coords_homogeneous[mask] @ T.T
        transformed_coords[mask, :] = transformed_coords_homogeneous[:, :-1]
        region_ids[mask] = entry.region_id

    adata.obs["x"] = transformed_coords[:, 0]
    adata.obs["y"] = transformed_coords[:, 1]
    adata.obs[REGION_FIELD] = region_ids
    return adata

def _get_transform_keys(e: "SpatialElement") -> list[str]:
    from spatialdata.transformations import get_transformation
    transformations = get_transformation(e, get_all=True)
    if not isinstance(transformations, dict):
        raise ValueError(f"This should be unreachable, get_transformation with get_all=True should always return a dict")
    return list(transformations.keys())

# ---- region resolution ----
def _resolve_regions_for_table(sdata: "SpatialData", table_name: str, sdata_name: str):
    from spatialdata.transformations import get_transformation, get_transformation_between_coordinate_systems
    from spatialdata.models import get_table_keys
    from anndata import AnnData
    adata = sdata.tables[table_name]
    if not isinstance(adata, AnnData):
        # if you hit this, it is likely because an invalid table_name was passed to this function.
        raise ValueError(f"Invalid table_name? '{table_name}' is not a table in '{sdata_name}':\n{sdata}")

    region, _obs_region_col, _instance_key = get_table_keys(adata)

    regions = region if isinstance(region, list) else [region]
    if not regions:
        raise ValueError(f"No region for table '{table_name}'")

    all_regions: dict[str, dict] = {}

    region_to_image: dict[str, ImageEntry] = {}
    
    for r in regions:
        annotated = sdata[r]  # shapes/points/labels/image ... probably shapes, but we need to find an appropriate image.
        if isinstance(annotated, AnnData):
            raise ValueError(f"This should be unreachable - '{r}' is not a table in '{sdata_name}':\n{sdata}")
        transformations = get_transformation(annotated, get_all=True)
        if not isinstance(transformations, dict):
            raise ValueError(f"This should be unreachable, get_transformation with get_all=True should always return a dict")
        cs_names = list(transformations.keys())
        if not cs_names:
            raise ValueError(f"Annotated element '{r}' has no coordinate system definitions")

        # Heuristic: prefer a CS literally named like the region, else the first
        cs = r if r in cs_names else cs_names[0]

        # images that are likely compatible
        img_candidates = [img for img in sdata.images.items() if cs in _get_transform_keys(img[1])]


        img_entries: list[ImageEntry] = []
        best_idx = 0
        for i, (img_path, img_obj) in enumerate(img_candidates):
            T = get_transformation_between_coordinate_systems(sdata, annotated, img_obj)
            if T is None:
                continue
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

            img_entries.append(ImageEntry(
                path=img_path,
                transform_to_image=T,
                extent_px=wh,
                is_primary=False,
                region_id=f"{sdata_name}_{r}"
            ))
            # primary selection heuristic: first viable with largest area
            prev = img_entries[best_idx].extent_px
            if wh and prev and (wh[0]*wh[1] > prev[0]*prev[1]):
                best_idx = len(img_entries)-1


        for j, entry in enumerate(img_entries):
            entry.is_primary = (j == best_idx)

        # todo - also add images that aren't the primary image...
        best_img = img_entries[best_idx]
        region_id = best_img.region_id
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
        
        region_to_image[r] = best_img

    # transform coordinates to image coordinates for each row in the table
    _transform_table_coordinates(adata, region_to_image)
    
    # store in uns for this table
    adata.uns.setdefault("mdv", {})
    adata.uns["mdv"].setdefault("regions", all_regions)


@dataclass
class SpatialDataConversionArgs:
    spatialdata_path: str
    output_folder: str
    preserve_existing: bool = False
    serve: bool = False

def _try_read_zarr(path: str):# -> sd.SpatialData | None:
    try:
        import spatialdata as sd
        return sd.read_zarr(path)
    except Exception as e:
        print(f"Warning: Failed to read SpatialData object from {path}: '{e}'")
        return None

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

    # from mdvtools.spatial.spatial_conversion import convert_spatialdata_to_mdv
    from mdvtools.conversions import convert_scanpy_to_mdv

    # we could do a nicer glob thing here, but I don't want to test that right now.
    sdata_paths = [
        os.path.join(args.spatialdata_path, f)
        for f in os.listdir(args.spatialdata_path)
    ]
    # sdata_paths = [f for f in sdata_paths if f.endswith(".zarr")]
    sdata_paths = sorted(sdata_paths)
    assert len(sdata_paths) > 0, "No SpatialData objects found in the folder"
    sdata_objects: dict[str, SpatialData] = {}
    adata_objects: list[AnnData] = []
    all_regions: dict[str, dict] = {}
    names: set[str] = set()
    for sdata_path in sdata_paths:
        sdata_name = os.path.basename(sdata_path)
        if sdata_name in names:
            raise ValueError(
                f"SpatialData object '{sdata_path}' has the same name as another object - this is not yet supported."
            )
        sdata = _try_read_zarr(sdata_path)
        if sdata is None:
            continue
        names.add(sdata_name)
        sdata_objects[sdata_name] = sdata
        # if there are no tables, maybe we could still try to get images as regions, just without any rows corresponding to them.
        for table_name, adata in sdata.tables.items():
            # this will have all kinds of side-effects on the adata
            # once it returns, it will have
            # - transformed coordinates in obs[x, y]
            # - regions metadata in uns["mdv"]["regions"]
            _resolve_regions_for_table(sdata, table_name, sdata_name)
            adata_objects.append(adata)
            # should these go in uns rather than obs?
            adata.obs["spatialdata_path"] = sdata_name
            adata.obs["table_name"] = table_name

            all_regions.update(adata.uns["mdv"]["regions"])

    if len(adata_objects) == 0:
        raise ValueError("No tables found in any SpatialData objects - this is not yet supported.")
    
    ## todo - try to make sure we have sparse CSC matrices if possible.
    merged_adata = ad_concat(adata_objects, index_unique="_")
    mdv = convert_scanpy_to_mdv(
        args.output_folder, merged_adata, delete_existing=not args.preserve_existing
    )
    os.makedirs(
        f"{mdv.dir}/spatial", exist_ok=True
    )  # pretty sure sdata.write will do this anyway
    for sdata_path, sdata in sdata_objects.items():
        sdata_name = os.path.basename(sdata_path)
        sdata.write(os.path.join(mdv.dir, "spatial", sdata_name))
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

    if args.serve:
        print(f"Serving project at {args.output_folder}")
        mdv.serve()
    else:
        print(f"Project saved to {args.output_folder}")

    return mdv
if __name__ == "__main__":
    # take an array of spatialdata objects from a folder and convert them to mdv.
    parser = argparse.ArgumentParser(description="Convert SpatialData to MDV format")
    parser.add_argument("spatialdata_path", type=str, help="Path to SpatialData data")
    parser.add_argument("output_folder", type=str, help="Output folder for MDV project")
    parser.add_argument("--preserve-existing", action="store_false", default=False, help="Preserve existing project data")
    parser.add_argument("--serve", action="store_false", default=False, help="Serve the project after conversion")
    args = parser.parse_args()

    mdv = convert_spatialdata_to_mdv(args) # type: ignore argparse->SpatialDataConversionArgs
