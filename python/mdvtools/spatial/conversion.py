import argparse
import os
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    import spatialdata as sd

def _find_default_image_for_cs(sdata: sd.SpatialData, coord_system="global"):
    sd = sdata.filter_by_coordinate_system(coord_system)
    if not hasattr(sd, "images"):
        return None
    items = list(sd.images.items())  
    return items[0][0] if items else None


if __name__ == "__main__":
    # take an array of spatialdata objects from a folder and convert them to mdv.
    parser = argparse.ArgumentParser(description="Convert SpatialData to MDV format")
    parser.add_argument("spatialdata_path", type=str, help="Path to SpatialData data")
    parser.add_argument("output_folder", type=str, help="Output folder for MDV project")
    parser.add_argument("--preserve-existing", action="store_false", default=False, help="Preserve existing project data")
    parser.add_argument("--serve", action="store_false", default=False, help="Serve the project after conversion")
    args = parser.parse_args()

    # imports can be slow, so doing them here rather than at the top of the file
    import anndata as ad
    import numpy as np
    import spatialdata as sd
    from spatialdata.transformations import get_transformation
    from spatialdata.models import get_table_keys
    # from mdvtools.spatial.spatial_conversion import convert_spatialdata_to_mdv
    from mdvtools.conversions import convert_scanpy_to_mdv

    sdata_paths = [os.path.join(args.spatialdata_path, f) for f in os.listdir(args.spatialdata_path)]
    sdata_paths = [f for f in sdata_paths if f.endswith(".zarr")]
    sdata_paths = sorted(sdata_paths)
    assert len(sdata_paths) > 0, "No SpatialData objects found in the folder"
    sdata_objects: dict[str, sd.SpatialData] = {}
    adata_objects: dict[str, ad.AnnData] = {}
    all_regions: dict[str, dict] = {}
    names: set[str] = set()
    for sdata_path in sdata_paths:
        sdata_name = os.path.basename(sdata_path)
        if sdata_name in names:
            raise ValueError(f"SpatialData object '{sdata_path}' has the same name as another object - this is not yet supported.")
        names.add(sdata_name)
        sdata = sd.read_zarr(sdata_path)
        sdata_objects[sdata_name] = sdata
        # FOR NOW::: assert that they each have a single coordinate system, image & table.
        if "table" not in sdata.tables:
            # pending implementation of support for other tables etc.
            raise ValueError(f"No default table found in SpatialData object '{sdata_path}' - this is not yet supported.")
        adata = sdata.tables["table"]
        adata.obs["spatialdata_path"] = sdata_name
        adata_objects[sdata_path] = adata
        # get the transformation for the image
        main_image = None
        main_image_name = None
        for image_name, image in sdata.images.items():
            main_image = image
            main_image_name = image_name
            print(f"Found image in sdata: {image_name}")
            break
        if len(sdata.images.items()) > 1:
            print(f"Warning: Multiple images found in SpatialData object '{sdata_path}' - this is not yet properly supported and may result in unexpected behaviour.")
        if main_image is None:
            # it should be valid to have no image - and w
            raise ValueError(f"No image found in SpatialData object '{sdata_path}' - this is not yet supported.")
        region, _element_description, _instance_key = get_table_keys(adata)
        transformation = get_transformation(sdata[region])
        if transformation is None:
            print(f"Warning: No transformation found for region {region} in SpatialData object '{sdata_path}' - this is unexpected, using Identity.")
            transformation = sd.transformations.Identity()
        ## Apply transformation matrix to spatial coordinates
        # Convert spatialdata transformation to affine matrix and apply to coordinates
        try:
            if isinstance(transformation, sd.transformations.Identity):
                # no point doing a load of matrix multiplication if the transformation is Identity
                coords = adata.obsm["spatial"]
                adata.obs["x"] = coords[:, 0]
                adata.obs["y"] = coords[:, 1]
            else:            
                # Get the affine matrix from the transformation
                # For spatial coordinates, we typically have x, y axes
                input_axes = ["x", "y"]
                output_axes = ["x", "y"]
                affine_matrix = transformation.to_affine_matrix(input_axes=input_axes, output_axes=output_axes)
                # nb - in the case of xenium, the transormation is Identity but we know there is a scale factor...
                
                # Convert coordinates to homogeneous coordinates (add 1s for translation)
                coords_homogeneous = np.column_stack([adata.obsm["spatial"], np.ones(adata.obsm["spatial"].shape[0])])            
                # Apply transformation: each row of coords_homogeneous @ affine_matrix
                transformed_coords_homogeneous = coords_homogeneous @ affine_matrix.T
                # Convert back from homogeneous coordinates (remove the last column)
                coords = transformed_coords_homogeneous[:, :-1]
                
                # we're keeping the original spatial coordinates in the obsm and adding these as x,y in obs
                # which will be used for position fields in the cells region data
                # in future we'll be able to use the untransformed coordinates, with the transform applied in shader
                # - we might be able to have a column data-type for 2d/3d coordinates and reduce marshalling in gpu buffer creation.
                adata.obs["x"] = coords[:, 0]
                adata.obs["y"] = coords[:, 1]
            # add a column to the obs to indicate the coordinate system
            # when we handle multiple coordinate systems within an sdata this will need to be updated.
            adata.obs["coordinate_system"] = sdata_name
            all_regions[sdata_name] = {
                # "width": main_image.shape[1],
                # "height": main_image.shape[0],
                "spatial": {
                    "coordinate_system": "global",
                    "file": sdata_name
                },
                "viv_image": {
                    # relative the "avivator" base set later
                    "file": f"{sdata_name}/images/{main_image_name}"
                },
                #!!! stuff we don't really use but liable to get a runtime error if we don't have.
                "roi": {
                    "min_x": 0,
                    "min_y": 0,
                    "max_x": 1000,
                    "max_y": 1000,
                },
                "images": {},
            }
        except Exception as e:
            print(f"Warning: Failed to transform spatial coordinates for {sdata_path}: '{e}'")
        # adata.obsm["spatial"] = coords

    # merged_adata = ad.concat(adata_objects.values(), label="coordinate_system", index_unique="_")
    merged_adata = ad.concat(adata_objects.values(), index_unique="_")
    mdv = convert_scanpy_to_mdv(args.output_folder, merged_adata, delete_existing=not args.preserve_existing)
    os.makedirs(f"{mdv.dir}/spatial", exist_ok=True) # pretty sure sdata.write will do this anyway
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
        "position_fields": ["x", "y"],
        "region_field": "coordinate_system",
        "default_color": "x", # todo: figure out a better way to determine this.
        "scale_unit": "µm",
        "scale": 1.0,
        "avivator": {
            "default_channels": [],
            "base_url": "spatial/",
        },
        "all_regions": all_regions
    }
    mdv.set_datasource_metadata(cells_md)


    if args.serve:
        print(f"Serving project at {args.output_folder}")
        mdv.serve()
    else:
        print(f"Project saved to {args.output_folder}")
