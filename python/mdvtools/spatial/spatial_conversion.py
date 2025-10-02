import argparse
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from spatialdata.models import SpatialElement
from mdvtools.spatial.xenium import convert_xenium_to_mdv
import spatialdata as sd
from spatialdata.transformations import (
    get_transformation,
    Scale,
    Identity,
)
from spatialdata.models import get_table_keys
from mdvtools.mdvproject import MDVProject
from mdvtools.conversions import get_matrix, _add_dims

# for the current implementation, we will add a column to the spatial datasource
# to indicate which coordinate system each cell belongs to.
# Values in this field should be able to identify which spatialdata store & coordinate system each cell belongs to.
REGION_FIELD = "spatial_coordinate_system"

def _find_default_image_for_cs(sdata: sd.SpatialData, coord_system = "global"):
    sd = sdata.filter_by_coordinate_system(coord_system)
    if not hasattr(sd, "images"):
        return None
    return sd.images.items()[0][0] # type: ignore

def _add_spatial_cs_to_project(mdv: MDVProject, coord_system: str, sdata: sd.SpatialData, sdata_path: str, default_color: str = REGION_FIELD):
    """
    Add a coordinate system as a region to the MDV project metadata.
    """

    ds = mdv.get_datasource_metadata("spatial")
    if ds is None:
        raise ValueError("No spatial datasource found")
    if 'regions' not in ds:
        ds['regions'] = {
            'region_field': REGION_FIELD, # this is a strong convention and mostly for internal use.
            'default_color': default_color, # in a script which is bespoke for a given dataset, we should be able to override this.
            # (for example, if we know that there's a 'cell_type' column in the table, we could set this to 'cell_type')
            
            'position_fields': ['x', 'y'], # ... and z, etc.
            
            # these shouldn't be specified at this level of the metadata - *we should deprecate these*.
            'scale_unit': 'µm', 
            'scale': 1.0,

    
            'avivator': {
                'default_channels': [],
                'base_url': 'spatial', # if we want remote data we need to review this
            },

            'all_regions': {}
        }
    

    region_name = f"{sdata_path}_{coord_system}" # if there is only one coordinate system, we should use the path as the region name.
    # we deliberately defer most processing to hypothetical future spatialdata.js runtime,
    # we just want to record that there is some spatial data to be found in the project.

    all_regions = ds['regions']['all_regions']

    image = _find_default_image_for_cs(sdata, coord_system)
    if image is None:
        raise ValueError(f"No image found for coordinate system {coord_system}")
    image_path = f"{sdata_path}/images/{image}"
    
    all_regions[region_name] = {
        "viv_image": {
            "file": image_path,
        },
        "spatial": {
            "coordinate_system": coord_system,
            "file": sdata_path,
        }
    }


def convert_multiple_spatialdata_to_mdv(
    folder: str,
    sdata_paths: list[str],
    max_dims: int = 3,
    delete_existing: bool = False,
    link_data: bool = False,
) -> MDVProject:
    
    mdv = MDVProject(folder, delete_existing=delete_existing)
    mdv.set_editable(True)

    for sdata_path in sdata_paths:
        sdata = sd.read_zarr(sdata_path)
        sdata_key = sdata_path.split('/')[-1]
        local_path = f"spatial/{sdata_key}"
        # can we do this in the background while we get on with the rest of the work?
        if not link_data: sdata.write(f"{mdv.dir}/{local_path}")
        else:
            # make a symlink to the data in the project/spatial/ directory.
            ...
        for coord_system in sdata.coordinate_systems:
            # how much work should we do in this function?
            # where do we merge our datasources?
            # where do we add our REGION_FIELD?
            _add_spatial_cs_to_project(mdv, coord_system, sdata, local_path)


    return mdv


def convert_spatialdata_to_mdv(
    folder: str,
    sdata: sd.SpatialData | str,
    max_dims: int = 3,
    delete_existing: bool = False,
    chunk_data: bool = True,
) -> MDVProject:
    """
    Convert a SpatialData object to an MDV project.
    
    Assumes that the SpatialData object is a single-cell dataset, with an AnnData table used to create cells & genes.
    
    TODO: column grouping, adding extra metadata, etc etc.
    """

    if isinstance(sdata, str):
        sdata = sd.read_zarr(sdata)
    mdv = MDVProject(folder, delete_existing=delete_existing)
    mdv.set_editable(True)

    adata = sdata.tables["table"]
    assert adata is not None # can we make this assumption?
    if adata.n_obs == 0 or adata.n_vars == 0:
        raise ValueError("Cannot convert empty AnnData (0 cells or 0 genes)")
    
    if len(sdata.coordinate_systems) > 1:
        # in future we will have nice ways of relating elements & coordinate systems
        # where possible, this should be left to spatialdata.js runtime rather than baked into mdv (meta)data.
        print(f"Warning: SpatialData has {len(sdata.coordinate_systems)} coordinate systems, we need to handle this better.")
    
    """
    The first element returned gives information regarding which spatial elements are annotated by the table, the second
    element gives information which column in table.obs contains the information which spatial element is annotated
    by each row in the table and the instance key indicates the column in obs giving information of the id of each row.
    """
    regions, region_key, instance_key = get_table_keys(adata)
    # e.g. ('202307141249_JaninaRun5-14-7-23_VMSC06901_Output_region_1_polygons', 'cells_region', 'EntityID')
    regions = [regions] if isinstance(regions, str) else regions
    cells_transform = Identity() # !!!hacky coordinate system handling
    for region in regions:
        element = sdata.get(region)
        assert isinstance(element, SpatialElement), f"Expected get_table_keys to always return regions relating to SpatialElement, got {type(element)} for {region}"
        t = get_transformation(element)
        if isinstance(t, dict):
            # todo: handle this
            print(f"Warning: Transform for '{region}' is a dict, which is not supported yet.")
            continue
        # this is bad, but if we just have one region we might get away with using this for cell coordinates.
        cells_transform = t
    
        
        
    
    print("Creating cells datasource")
    cell_table = adata.obs.copy()
    cell_table["cell_id"] = cell_table.index
    
    if "spatial" in adata.obsm:
        spatial_coords = adata.obsm["spatial"].copy()
        
        region = regions[0]
        element = sdata.get(region)
        assert isinstance(element, SpatialElement), f"Expected get_table_keys to always return regions relating to SpatialElement, got {type(element)} for {region}"
        t = get_transformation(element)
        if isinstance(t, dict):
            # todo: handle this
            raise ValueError(f"Transform for '{region}' is a dict, which is not supported yet.")
            
        
        # rather than have a single cells_transform, each row should be transformed corresponding to the adata.obs[region_key]
        # ultimately we probably don't want to bake these kinds of transforms into the mdv data... 
        # but ultimately, all the data should be read from the spatialdata store at runtime, along with appropriate transforms.
        spatial_coords = sd.transform(spatial_coords, cells_transform)
        
        cell_table["x"] = spatial_coords[:, 0]
        cell_table["y"] = spatial_coords[:, 1]
        # we should also handle time when relevant, and make sure we have 3d/4d regions etc etc... not much good will come of the current approach.
        if spatial_coords.shape[1] == 3:
            cell_table["z"] = spatial_coords[:, 2]
    else:
        print("Warning: No spatial coordinates found in obsm['spatial']")
        # Try to find coordinates in other common locations

    # Add dimensionality reductions
    cell_table = _add_dims(cell_table, adata.obsm, max_dims)
    # todo: consider what happens if we are merging multiple datasets...
    # columns = [{"name": "cell_id", "datatype": "unique"}]
    mdv.add_datasource(f"spatial", cell_table)

    # todo: something different.......
    label = ""

    # Create genes datasource - nb, mock spatialdata has a varm that doesn't really represent genes (channel_0_sum...)
    print("Creating genes datasource")
    gene_table = adata.var.copy()
    gene_table["gene_id"] = gene_table.index
    gene_table = _add_dims(gene_table, adata.varm, max_dims)
    mdv.add_datasource(f"{label}genes", gene_table)

    # Link cells and genes through expression data
    print("Linking cells and genes")
    mdv.add_rows_as_columns_link(f"spatial", f"{label}genes", "gene_id", "Gene Expr")

    # Add gene expression matrix
    print("Adding gene expression matrix")
    # todo: add another column for density (so that we can identify genes that are less sparse in the data)?
    matrix, sparse = get_matrix(adata.X)
    if matrix is not None and matrix.shape[1] != 0: # type: ignore wtf makes this necessary?
        mdv.add_rows_as_columns_subgroup(
            f"spatial", f"{label}genes", "gs", matrix, 
            name="gene_scores", label="Gene Scores",
            chunk_data=chunk_data
        )
    
    # Add additional layers if present
    for layer_name, layer_matrix in adata.layers.items():
        print(f"Adding layer {layer_name}")
        layer_matrix, layer_sparse = get_matrix(layer_matrix)
        if layer_matrix is not None and layer_matrix.shape[1] != 0: # type: ignore
            mdv.add_rows_as_columns_subgroup(
                f"spatial", f"{label}genes", layer_name, layer_matrix,
                chunk_data=chunk_data
            )
    
    # Set up spatial region data for cells
    for cs in sdata.coordinate_systems:
        # it isn't clear that this is what anyone would want or expect, but I think the best mapping to mdv as of writing
        # is that each "coordinate system" is a region, with any related SpatialElements belonging to that.
        # anticipating a major revision of how spatialdata gets represented, 
        # which may mean our current 'regions' metadata is not relevant.    
        # Previously, I'd been tending to use each "region" as an image, rather than an ROI within it.
        # in spatialdata, it seems that if you have e.g. an image, and shapes for cells etc, 
        # the coordinate systems is the thing that relates things pertaining to that image/sample.
        ...
        
    # if "x" in cell_table.columns and "y" in cell_table.columns:
    #     print("Setting up spatial region data")
    #     try:
    #         mdv.set_region_data(
    #             f"spatial",
    #             cell_table,
    #             region_field="cell_id",
    #             default_color="total_counts" if "total_counts" in cell_table.columns else "cell_id",
    #             position_fields=["x", "y"],
    #             scale_unit="µm",
    #             scale=1.0 # todo: adjust based on provided scale info. in future, more sopbisticated coordinate system handling
    #         )
    #     except Exception as e:
    #         print(f"Warning: Could not set up region data: {str(e)}")
    

    # todo: current idea is that these won't be DataSources, but things that are associated with cells regions
    #  and understood by SpatialLayers
    # we include an entire spatialdata.zarr store - and aim to use this for as much as possible in future...
    # a project should be able to reference multiple spatialdata stores either local to the project, remote, 
    # or in some other local place that may be shared with other projects.
    sdata.write(f"{mdv.dir}/spatial/spatialdata.zarr")
    return mdv



# if __name__ == "__main__X":
#     sdata = sd.datasets.blobs()
#     mdv = convert_spatialdata_to_mdv(
#         folder="/app/mdv/mock_sdata_blobs",
#         spatialdata_path="/app/mdv/mock_sdata_blobs/spatialdata.zarr",
#         delete_existing=True,
#     )
#     mdv.serve()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert SpatialData to MDV format")
    parser.add_argument("spatialdata_path", type=str, help="Path to SpatialData data")
    parser.add_argument("output_folder", type=str, help="Output folder for MDV project")
    parser.add_argument("--delete-existing", action="store_false", default=True, help="Delete existing project data")
    args = parser.parse_args()

    mdv = convert_spatialdata_to_mdv(
        folder=args.output_folder,
        sdata=args.spatialdata_path,
        delete_existing=args.delete_existing,
    )
    mdv.serve()
