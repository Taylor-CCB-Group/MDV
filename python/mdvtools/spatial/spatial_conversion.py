import argparse
from typing import Optional, List
from mdvtools.spatial.xenium import convert_xenium_to_mdv
import spatialdata as sd
from spatialdata.transformations import (
    get_transformation,
    Scale,
)
from spatialdata_io import xenium
import copy
import subprocess
import os
import tempfile
import shutil

from mdvtools.mdvproject import MDVProject
from mdvtools.conversions import get_matrix, _add_dims

def convert_spatialdata_to_mdv(
    folder: str,
    sdata: sd.SpatialData | str,
    max_dims: int = 3,
    delete_existing: bool = False,
    label: str = "",
    chunk_data: bool = False,
) -> MDVProject:
    if isinstance(sdata, str):
        sdata = sd.read_zarr(sdata)
    mdv = MDVProject(folder, delete_existing=delete_existing)
    mdv.set_editable(True)

    adata = sdata.tables["table"]
    assert adata is not None # can we make this assumption?
    if adata.n_obs == 0 or adata.n_vars == 0:
        raise ValueError("Cannot convert empty AnnData (0 cells or 0 genes)")
    
    print("Creating cells datasource")
    cell_table = adata.obs.copy()
    cell_table["cell_id"] = cell_table.index
    if "spatial" in adata.obsm:
        spatial_coords = adata.obsm["spatial"]
        cell_table["x"] = spatial_coords[:, 0]
        cell_table["y"] = spatial_coords[:, 1]
        # consider non-2d here...
    else:
        print("Warning: No spatial coordinates found in obsm['spatial']")
        # Try to find coordinates in other common locations

    # Add dimensionality reductions
    cell_table = _add_dims(cell_table, adata.obsm, max_dims)
    # todo: consider what happens if we are merging multiple datasets...
    # columns = [{"name": "cell_id", "datatype": "unique"}]
    mdv.add_datasource(f"{label}cells", cell_table)

    # Create genes datasource - nb, mock spatialdata has a varm that doesn't really represent genes (channel_0_sum...)
    print("Creating genes datasource")
    gene_table = adata.var.copy()
    gene_table["gene_id"] = gene_table.index
    gene_table = _add_dims(gene_table, adata.varm, max_dims)
    mdv.add_datasource(f"{label}genes", gene_table)

    # Link cells and genes through expression data
    print("Linking cells and genes")
    mdv.add_rows_as_columns_link(f"{label}cells", f"{label}genes", "gene_id", "Gene Expr")

    # Add gene expression matrix
    print("Adding gene expression matrix")
    matrix, sparse = get_matrix(adata.X)
    if matrix is not None and matrix.shape[1] != 0: # type: ignore wtf makes this necessary?
        mdv.add_rows_as_columns_subgroup(
            f"{label}cells", f"{label}genes", "gs", matrix, 
            name="gene_scores", label="Gene Scores",
            chunk_data=chunk_data
        )
    
    # Add additional layers if present
    for layer_name, layer_matrix in adata.layers.items():
        print(f"Adding layer {layer_name}")
        layer_matrix, layer_sparse = get_matrix(layer_matrix)
        if layer_matrix is not None and layer_matrix.shape[1] != 0: # type: ignore
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
    
    # todo: handle images and other spatialdata properly

    # todo: current idea is that these won't be DataSources, but things that are associated with cells regions
    #  and understood by SpatialLayers
    # we include an entire spatialdata.zarr store - and aim to use this for as much as possible in future...
    sdata.write(f"{mdv.dir}/spatialdata.zarr")
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
