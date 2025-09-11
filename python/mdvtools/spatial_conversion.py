from typing import Optional, List
import spatialdata as sd
from spatialdata_io import xenium
import copy
import subprocess
import os
import tempfile
import shutil

from mdvtools.mdvproject import MDVProject
from mdvtools.conversions import get_matrix, _add_dims


def convert_xenium_to_mdv(
    folder: str,
    xenium_path: str,
    max_dims: int = 3,
    delete_existing: bool = False,
    label: str = "",
    chunk_data: bool = False,
    # todo handling of cell boundaries, transcripts, transcripts.
    load_cell_boundaries: bool = False,
    load_nucleus_boundaries: bool = False,
    load_transcripts: bool = False,
    load_morphology_images: bool = True,
    load_cells_table: bool = True,
    selected_genes: Optional[List[str]] = None,
    gene_identifier_column: Optional[str] = None,
    n_jobs: int = 1
) -> MDVProject:
    """
    Convert Xenium spatial transcriptomics data to MDV (Multi-Dimensional Viewer) format.

    This function uses spatialdata-io to read Xenium data and transforms it into 
    the MDV project structure, handling cells, genes, spatial coordinates, and 
    optional additional data like cell boundaries, transcripts, and images.

    Args:
        folder (str): Path to the target MDV project folder
        xenium_path (str): Path to the Xenium dataset directory
        max_dims (int, optional): Maximum number of dimensions to include from 
            dimensionality reductions. Defaults to 3.
        delete_existing (bool, optional): Whether to delete existing project data. 
            If False, merges with existing data. Defaults to False.
        label (str, optional): Prefix to add to datasource names and metadata columns
            when merging with existing data. Defaults to "".
        chunk_data (bool, optional): For dense matrices, transposing and flattening
            will be performed in chunks. Saves memory but takes longer. Default is False.
        load_cell_boundaries (bool, optional): Whether to load cell boundary polygons. 
            Defaults to True.
        load_nucleus_boundaries (bool, optional): Whether to load nucleus boundary polygons. 
            Defaults to True.
        load_transcripts (bool, optional): Whether to load individual transcript locations. 
            Can be memory intensive for large datasets. Defaults to False.
        load_morphology_images (bool, optional): Whether to load morphology images. 
            Defaults to True.
        load_cells_table (bool, optional): Whether to load the cell annotations table. 
            Defaults to True.
        selected_genes (List[str], optional): List of specific genes to load. 
            If None, loads all genes. Defaults to None.
        gene_identifier_column (str, optional): Column name for gene identifiers. 
            If None, uses 'name' column. Defaults to None.
        n_jobs (int, optional): Number of parallel jobs for reading data. Defaults to 1.

    Returns:
        MDVProject: The configured MDV project object with the converted data

    Notes:
        Data Structure Creation:
        - Creates main datasources: '{label}cells' and '{label}genes'
        - Preserves all cell metadata from Xenium cell table
        - Preserves all gene metadata from Xenium features
        - Adds spatial coordinates as x, y columns
        - Links cells and genes through expression data
        - Optionally adds cell/nucleus boundaries as separate datasources
        - Optionally adds transcript locations as points datasource
        - Optionally adds morphology images as OME-NGFF image datasources

        Spatial Data Handling:
        - Sets up region data for spatial visualization
        - Configures position fields for spatial plotting
        - Handles coordinate system transformations

    Raises:
        ValueError: If the Xenium data cannot be loaded or is invalid
        IOError: If there are issues with file operations
        Exception: For other unexpected errors during conversion
    """
    
    # Load Xenium data using spatialdata-io
    print(f"Loading Xenium data from {xenium_path}")
    try:
        sdata: sd.SpatialData = xenium(
            xenium_path,
            cells_boundaries=load_cell_boundaries,
            nucleus_boundaries=load_nucleus_boundaries,
            transcripts=load_transcripts,
            morphology_focus=load_morphology_images,
            cells_table=load_cells_table,
            n_jobs=n_jobs
        )
    except Exception as e:
        raise ValueError(f"Failed to load Xenium data: {str(e)}")
    
    # Validate that we have the required data
    if "table" not in sdata.tables:
        raise ValueError("No cell table found in Xenium data")
    
    adata = sdata.tables["table"]
    if adata.n_obs == 0 or adata.n_vars == 0:
        raise ValueError("Cannot convert empty Xenium data (0 cells or 0 genes)")
    
    # Filter genes if specified
    if selected_genes is not None:
        # Find genes that exist in the dataset
        available_genes = adata.var_names.intersection(selected_genes)
        if len(available_genes) == 0:
            raise ValueError("None of the selected genes found in the dataset")
        adata = adata[:, available_genes].copy()
        print(f"Filtered to {len(available_genes)} selected genes")
    
    # Initialize MDV project
    mdv = MDVProject(folder, delete_existing=delete_existing)
    
    # Preserve current views if not deleting existing
    current_views = None
    if not delete_existing:
        current_views = mdv.views
    else:
        current_views = {}
    
    # Create cells datasource
    print("Creating cells datasource")
    cell_table = adata.obs.copy()
    cell_table["cell_id"] = cell_table.index
    
    # Add spatial coordinates if available
    if "spatial" in adata.obsm:
        spatial_coords = adata.obsm["spatial"]
        cell_table["x"] = spatial_coords[:, 0]
        cell_table["y"] = spatial_coords[:, 1]
    else:
        print("Warning: No spatial coordinates found in obsm['spatial']")
        # Try to find coordinates in other common locations
        if hasattr(adata, 'obsm') and len(adata.obsm) > 0:
            for key, coords in adata.obsm.items():
                if coords.shape[1] >= 2:
                    cell_table["x"] = coords[:, 0]
                    cell_table["y"] = coords[:, 1]
                    print(f"Using coordinates from obsm['{key}']")
                    break
    
    # Add dimensionality reductions
    cell_table = _add_dims(cell_table, adata.obsm, max_dims)
    
    # Define cell_id as unique column
    columns = [{"name": "cell_id", "datatype": "unique"}]
    mdv.add_datasource(f"{label}cells", cell_table, columns)
    
    # Create genes datasource
    print("Creating genes datasource")
    gene_table = adata.var.copy()
    
    # Set up gene identifier column
    if gene_identifier_column and gene_identifier_column not in gene_table.columns:
        print(f"Gene identifier column {gene_identifier_column} not found, using index")
        gene_identifier_column = None
    if not gene_identifier_column:
        gene_identifier_column = f"{label}name"
        gene_table[gene_identifier_column] = gene_table.index
    
    # Add dimensionality reductions for genes if available
    gene_table = _add_dims(gene_table, adata.varm, max_dims)
    
    mdv.add_datasource(f"{label}genes", gene_table)
    
    # Link cells and genes through expression data
    print("Linking cells and genes")
    mdv.add_rows_as_columns_link(f"{label}cells", f"{label}genes", gene_identifier_column, "Gene Expr")
    
    # Add gene expression matrix
    print("Adding gene expression matrix")
    matrix, sparse = get_matrix(adata.X)
    if matrix is not None and matrix.shape[1] != 0:
        mdv.add_rows_as_columns_subgroup(
            f"{label}cells", f"{label}genes", "gs", matrix, 
            name="gene_scores", label="Gene Scores",
            chunk_data=chunk_data
        )
    
    # Add additional layers if present
    for layer_name, layer_matrix in adata.layers.items():
        print(f"Adding layer {layer_name}")
        layer_matrix, layer_sparse = get_matrix(layer_matrix)
        if layer_matrix is not None and layer_matrix.shape[1] != 0:
            mdv.add_rows_as_columns_subgroup(
                f"{label}cells", f"{label}genes", layer_name, layer_matrix,
                chunk_data=chunk_data
            )
    
    # Set up spatial region data for cells
    if "x" in cell_table.columns and "y" in cell_table.columns:
        print("Setting up spatial region data")
        try:
            mdv.set_region_data(
                f"{label}cells",
                cell_table,
                region_field="cell_id",
                default_color="total_counts" if "total_counts" in cell_table.columns else "cell_id",
                position_fields=["x", "y"],
                scale_unit="µm",
                scale=1.0 # todo: adjust based on provided scale info. in future, more sopbisticated coordinate system handling
            )
        except Exception as e:
            print(f"Warning: Could not set up region data: {str(e)}")
    
    # todo: current idea is that these won't be DataSources, but things that are associated with cells regions
    #  and understood by SpatialLayers
    # Add cell boundaries if available
    # if load_cell_boundaries and "cell_boundaries" in sdata.shapes:
    #     print("Adding cell boundaries")
    #     try:
    #         cell_boundaries = sdata.shapes["cell_boundaries"]
    #         # Convert GeoDataFrame to DataFrame for MDV
    #         boundaries_df = pd.DataFrame({
    #             "cell_id": cell_boundaries.index,
    #             "geometry": cell_boundaries.geometry
    #         })
    #         mdv.add_datasource(f"{label}cell_boundaries", boundaries_df)
    #     except Exception as e:
    #         print(f"Warning: Could not add cell boundaries: {str(e)}")
    
    # # Add nucleus boundaries if available
    # if load_nucleus_boundaries and "nucleus_boundaries" in sdata.shapes:
    #     print("Adding nucleus boundaries")
    #     try:
    #         nucleus_boundaries = sdata.shapes["nucleus_boundaries"]
    #         boundaries_df = pd.DataFrame({
    #             "nucleus_id": nucleus_boundaries.index,
    #             "geometry": nucleus_boundaries.geometry
    #         })
    #         mdv.add_datasource(f"{label}nucleus_boundaries", boundaries_df)
    #     except Exception as e:
    #         print(f"Warning: Could not add nucleus boundaries: {str(e)}")
    
    # # Add transcripts if available
    # if load_transcripts and "transcripts" in sdata.points:
    #     print("Adding transcripts")
    #     try:
    #         transcripts = sdata.points["transcripts"]
    #         # Convert to DataFrame
    #         transcripts_df = pd.DataFrame(transcripts)
    #         mdv.add_datasource(f"{label}transcripts", transcripts_df)
    #     except Exception as e:
    #         print(f"Warning: Could not add transcripts: {str(e)}")
    
    # Add morphology images as OME-NGFF format if available
    if load_morphology_images and "morphology_focus" in sdata.images:
        print("Adding morphology images as JP2K OME-TIFF")
        try:
            morphology_image = sdata.images["morphology_focus"]
            # _save_ome_ngff_image(mdv, morphology_image, f"{label}morphology_focus")
            _convert_morphology_image(mdv, xenium_path)
        except Exception as e:
            print(f"Warning: Could not add morphology images: {str(e)}")
    
    # Set up views
    if delete_existing:
        # Create new default view
        mdv.set_view("default", {"initialCharts": {"cells": [], "genes": []}}, True)
        mdv.set_editable(True)
    else:
        # Update existing views with new datasources
        new_views = {}
        for view_name, view_data in current_views.items():
            new_view_data = copy.deepcopy(view_data)
            
            # Initialize new charts if they don't exist
            if "initialCharts" not in new_view_data:
                new_view_data["initialCharts"] = {}
            
            # Add new datasources to initialCharts
            new_view_data["initialCharts"][f"{label}cells"] = []
            new_view_data["initialCharts"][f"{label}genes"] = []
            
            # Initialize dataSources if they don't exist
            if "dataSources" not in new_view_data:
                new_view_data["dataSources"] = {}
            
            # Add new datasources with panel widths
            new_view_data["dataSources"][f"{label}cells"] = {"panelWidth": 50}
            new_view_data["dataSources"][f"{label}genes"] = {"panelWidth": 50}
            
            new_views[view_name] = new_view_data
        
        mdv.views = new_views
    
    print("Xenium to MDV conversion completed successfully")
    return mdv


def convert_xenium_zarr_to_mdv(
    folder: str,
    zarr_path: str,
    max_dims: int = 3,
    delete_existing: bool = False,
    label: str = "",
    chunk_data: bool = False,
    selected_genes: Optional[List[str]] = None,
    gene_identifier_column: Optional[str] = None
) -> MDVProject:
    """
    Convert Xenium data from Zarr format to MDV format.
    
    This is a convenience function for when Xenium data has already been converted
    to Zarr format using spatialdata-io.
    
    Args:
        folder (str): Path to the target MDV project folder
        zarr_path (str): Path to the Zarr store containing Xenium data
        max_dims (int, optional): Maximum number of dimensions to include from 
            dimensionality reductions. Defaults to 3.
        delete_existing (bool, optional): Whether to delete existing project data. 
            If False, merges with existing data. Defaults to False.
        label (str, optional): Prefix to add to datasource names and metadata columns
            when merging with existing data. Defaults to "".
        chunk_data (bool, optional): For dense matrices, transposing and flattening
            will be performed in chunks. Saves memory but takes longer. Default is False.
        selected_genes (List[str], optional): List of specific genes to load. 
            If None, loads all genes. Defaults to None.
        gene_identifier_column (str, optional): Column name for gene identifiers. 
            If None, uses 'name' column. Defaults to None.

    Returns:
        MDVProject: The configured MDV project object with the converted data
    """
    
    print(f"Loading Xenium data from Zarr store: {zarr_path}")
    try:
        sdata = sd.read_zarr(zarr_path)
    except Exception as e:
        raise ValueError(f"Failed to load Zarr data: {str(e)}")
    
    # Extract the AnnData object
    if "table" not in sdata.tables:
        raise ValueError("No cell table found in Zarr data")
    
    adata = sdata.tables["table"]
    
    # Use the main conversion function with the loaded data
    return convert_xenium_to_mdv(
        folder=folder,
        xenium_path="",  # Not used since we already have the data
        max_dims=max_dims,
        delete_existing=delete_existing,
        label=label,
        chunk_data=chunk_data,
        load_cell_boundaries="cell_boundaries" in sdata.shapes,
        load_nucleus_boundaries="nucleus_boundaries" in sdata.shapes,
        load_transcripts="transcripts" in sdata.points,
        load_morphology_images="morphology_focus" in sdata.images,
        load_cells_table=True,
        selected_genes=selected_genes,
        gene_identifier_column=gene_identifier_column
    )

def _convert_morphology_image(mdv: MDVProject, xenium_path: str) -> None:
    # use bioformats2raw & raw2ometiff to convert the image to JP2K OME-TIFF
    # this is particularly flakey because we don't have those as dependencies - also require java, blosc, ...
    # first test is running locally with bespoke environment implied, may update docker container to include them, or not.
    # get a directory listing of xenium_path/morphology_focus, then find the one of the ome.tif files for bioformats2raw
    # then use raw2ometiff to convert the raw file to OME-TIFF
    file = os.listdir(f"{xenium_path}/morphology_focus")[0] # todo: more robust way to find the file
    with tempfile.TemporaryDirectory() as temp_dir:
        tmp_zarr = f"{temp_dir}/{file}.zarr"
        tmp_tiff = f"{temp_dir}/jp2_{file}"
        print(f"Converting {file} to JP2K OME-TIFF")
        print(f"tmp_zarr: '{tmp_zarr}'")
        print(f"tmp_tiff: '{tmp_tiff}'")
        subprocess.run(["bioformats2raw", f"{xenium_path}/morphology_focus/{file}", tmp_zarr])
        subprocess.run([
            "raw2ometiff", tmp_zarr, tmp_tiff, 
            "--compression", "JPEG-2000 Lossy", "--quality", "10"
        ])
        # move the tiff file to the images directory
        os.makedirs(f"{mdv.dir}/images/avivator", exist_ok=True)
        shutil.move(tmp_tiff, f"{mdv.dir}/images/avivator/jp2_{file}")
        ds = mdv.get_datasource_metadata("cells")
        # First-pass, CBA to try to use the methods documented in spatialdata.md - I think we may want a refactor
        # this is somewhat based on mdv.update_datasource_for_tiff - also sus.
        # mdv.set_region_data("cells", region_field="region", default_color="x")
        # mdv.add_viv_images("cells", [])
        ds["regions"] = {
            "position_fields": ["x", "y"],
            "region_field": "region",
            "default_color": "region",
            "scale_unit": "µm",
            "scale": 1.0, # todo: read appropriate scale/roi from the image metadata
            "all_regions": {
                # definitely not wanting this as region id but looks like it might work for very initial testing
                # ultimately want to be able to have multiple e.g. xenium inputs in the same project
                # - so each input can have a corresponding region id... and we might have ways of loading column data from a given store...
                "cell_circles": {
                    "roi": {
                        "min_x": 0,
                        "min_y": 0,
                        "max_x": 100,
                        "max_y": 100
                    },
                    "images": {},
                    "viv_image": {
                        "file": f"jp2_{file}",
                        "linked_file": False,
                    }
                }
            },
            "avivator": {
                "default_channels": [],
                "base_url": "images/avivator"
            }
        }
        print(ds)
        mdv.set_datasource_metadata(ds)

if __name__ == "__main__":
    import os
    folder = os.path.expanduser("~/data/mdv_xenium_test")
    x_path = os.path.expanduser("~/data/Xenium_Prime_MultiCellSeg_Mouse_Ileum_tiny_outs")
    mdv = convert_xenium_to_mdv(
        folder=folder,
        xenium_path=x_path,
        delete_existing=True,
    )
    mdv.serve()