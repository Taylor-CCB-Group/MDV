from scanpy import AnnData
import scipy
import pandas as pd
from os.path import join, split
import os
import re
import json
import gzip
import copy
import yaml
import shutil
import h5py
import logging
from pathlib import Path
from typing import Dict, List, Tuple
from werkzeug.utils import secure_filename
import numpy as np
from .mdvproject import MDVProject,create_bed_gz_file

logger = logging.getLogger(__name__)


def convert_scanpy_to_mdv(
    folder: str, 
    scanpy_object: AnnData, 
    max_dims: int = 3, 
    delete_existing: bool = False, 
    label: str = "",
    chunk_data: bool = False,
    add_layer_data = True,
    gene_identifier_column = None
) -> MDVProject:
    """
    Convert a Scanpy AnnData object to MDV (Multi-Dimensional Viewer) format.

    This function transforms single-cell RNA sequencing data from AnnData format into 
    the MDV project structure, handling both cells and genes as separate datasources 
    with their associated dimensionality reductions and metadata.

    Args:
        folder (str): Path to the target MDV project folder
        scanpy_object (AnnData): The AnnData object containing the single-cell data
        max_dims (int, optional): Maximum number of dimensions to include from 
            dimensionality reductions. Defaults to 3.
        delete_existing (bool, optional): Whether to delete existing project data. 
            If False, merges with existing data. Defaults to False.
        label (str, optional): Prefix to add to datasource names and metadata columns
            when merging with existing data. Defaults to "".
        chunk_data (bool, optional): For dense matrices, transposing and flattening
            will be performed in chunks. Saves memory but takes longer. Default is False.
        add_layer_data (bool, optional): If True (default) then the layer data (log values etc.)
            will be added, otherwise just the X object will be used
        gene_identifier_column: (str, optional) This is the gene column that the user will use to
            identify the gene. If not specified (default) than a column 'name' will be added that is
            created from the index (which is usaully the unique gene name)
    Returns:
        MDVProject: The configured MDV project object with the converted data

    Notes:
        Data Structure Creation:
        - Creates two main datasources: '{label}cells' and '{label}genes'
        - Preserves all cell metadata from scanpy_object.obs
        - Preserves all gene metadata from scanpy_object.var
        - Transfers dimension reductions from obsm/varm matrices
        - Links cells and genes through expression data
        - Adds gene expression scores as a subgroup

        View Handling:
        - If delete_existing=True:
            * Creates new default view with empty initial charts
            * Sets project as editable
        - If delete_existing=False:
            * Preserves existing views
            * Updates views with new datasources
            * Maintains panel widths and other view settings
            * Adds new datasources to each view's initialCharts

        Dimension Reduction:
        - Processes dimensionality reductions up to max_dims
        - Supports standard formats (e.g., PCA, UMAP, t-SNE)
        - Column names in the format: {reduction_name}_{dimension_number}

    Raises:
        ValueError: If the provided AnnData object is invalid or missing required components
        IOError: If there are issues with file operations in the target folder
        Exception: For other unexpected errors during conversion
    """
    # Validate input AnnData
    if scanpy_object.n_obs == 0 or scanpy_object.n_vars == 0:
        raise ValueError("Cannot convert empty AnnData object (0 cells or 0 genes)")
    
    mdv = MDVProject(folder, delete_existing=delete_existing)

    # If not deleting existing, preserve current views
    current_views = None
    if not delete_existing:
        current_views = mdv.views
    else:
        current_views = {}
    # create datasource 'cells'
    cell_table = scanpy_object.obs
    cell_table["cell_id"] = cell_table.index

    # add any dimension reduction to the dataframe
    cell_table = _add_dims(cell_table, scanpy_object.obsm, max_dims)

    # cell_id is the unique barcode and should of type unique
    # (will be text16 by default if number of values are below 65536)
    # hopefully other columns will be of the correct format
    columns=[{
        "name":"cell_id",
        "datatype":"unique"
    }]
    mdv.add_datasource(f"{label}cells", cell_table,columns)

    # create datasource 'genes'
    gene_table = scanpy_object.var
    #need a way of detecting which column is the common gene name
    #most times it is the index but sometimes this is just the gene code or an
    #incremental number.
    if gene_identifier_column and not gene_identifier_column in gene_table.columns:
        print(f"gene identifier column {gene_identifier_column} not found, using index")
        gene_identifier_column= None
    if not gene_identifier_column:
        gene_identifier_column= f"{label}name"
        gene_table[gene_identifier_column] = gene_table.index
    gene_table = _add_dims(gene_table, scanpy_object.varm, max_dims)

    #originally column had to be unique - but now is just text
    #need to coerce gene name column to unique
    #columns=[{"name":"name","datatype":"text16"}]
    mdv.add_datasource(f"{label}genes", gene_table)

    # link the two datasets
    mdv.add_rows_as_columns_link(f"{label}cells", f"{label}genes", gene_identifier_column, "Gene Expr")

    #get the matrix in the correct format
    print("Getting Matrix")
    matrix,sparse= get_matrix(scanpy_object.X)
    #sometimes X is empty - all the data is in the layers
    assert matrix is not None # asserting here so that 'invalid_adata' test gets expected error
    # don't want to faff about with more complex logic for permutation with add_layer_data just now
    if matrix.shape[1] !=0: # type: ignore maybe change get_matrix in future so type inference is better
        # add the gene expression
        print("Adding gene expression")
        mdv.add_rows_as_columns_subgroup(
            f"{label}cells", f"{label}genes", "gs", matrix, name="gene_scores", label="Gene Scores",
            # sparse=sparse, #this should be inferred from the matrix
            chunk_data=chunk_data
        )

    #now add layers
    if add_layer_data:
        for layer,matrix in scanpy_object.layers.items():
            matrix,sparse = get_matrix(matrix)
            print(f"Adding layer {layer}")
            mdv.add_rows_as_columns_subgroup(
                f"{label}cells", f"{label}genes", layer, matrix, chunk_data=chunk_data
            )

    if delete_existing:
        # If we're deleting existing, create new default view
        mdv.set_view("default", {"initialCharts": {"cells": [], "genes": []}}, True)
        mdv.set_editable(True)
    else:
        # If we're not deleting existing, update existing views with new datasources
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
        
    return mdv

def convert_mudata_to_mdv(folder,mudata_object,max_dims=3,delete_existing=False, chunk_data=False):
    md=mudata_object
    p= MDVProject(folder,delete_existing=delete_existing)
    #are there any general drs
    table = _add_dims(md.obs,md.obsm,max_dims)
    #add drs derived from modalities
    for name,mod in md.mod.items():
        table = _add_dims(table,mod.obsm,max_dims,name)
    md.obs["cell_id"] = md.obs.index
    columns= [{"name":"cell_id","datatype":"unique"}]
    p.add_datasource("cells",table,columns)

    for mod in md.mod.keys():
        mdata = md.mod[mod]
        #add the modality to project
        p.add_datasource(mod,mdata.var)
        #adds the index to the data as a name column
        #This is usually the unique gene name - but not always
        column ="name"
        #no longer unique
        #column = {"name":"name","datatype":"unique"}
        p.set_column(mod,column,mdata.var.index)
        #mod is used as both the tag and the label
        #the name column is specified as the identifier that the user will use
        #it is derived from the index and is usually the gene 'name'
        #However it may not be appropriate can be changed later on 
        p.add_rows_as_columns_link("cells",mod,"name",mod)
        matrix,sparse= get_matrix(mdata.X)
        #sometimes X is empty - all the data is in the layers
        if matrix.shape[1] !=0:
            p.add_rows_as_columns_subgroup("cells",mod,mod+"_expr",matrix,sparse=sparse, chunk_data=chunk_data)
        #now add the layers (if there are any)
        layers = mdata.layers.keys()
        for layer in layers:
            matrix  = mdata.layers[layer]
            matrix,sparse = get_matrix(matrix,mdata.obs_names,md.obs_names)
            p.add_rows_as_columns_subgroup("cells",mod,f"{mod}_{layer}",matrix,sparse=sparse, chunk_data=chunk_data)
    return p


# The main_names correspond to obs_names in the main anndata object
# and the mod_names those in the modality.
# Only required if data is being added from a modality,as sometimes
# the modality's obs_names will be in a different order and/or a subset of the main names
# Hence a sparse matrix corresponding to the main indices needs to be created
def get_matrix(matrix,main_names=None,mod_names=None) -> tuple[scipy.sparse.csc_matrix | np.ndarray, bool]:
    if main_names is None:
        main_names = []
    if mod_names is None:
        mod_names = []
    #check where the matrix data actually is
    matrix = matrix.value if hasattr(matrix,"value") else matrix
    #is the matrix backed sparse matrix -convert to non backed else cannot convert
    if hasattr(matrix,"backend"):
        matrix=matrix._to_backed()
    #not a sparse matrix - nothing else to do
    if not scipy.sparse.issparse(matrix):
        return matrix,False
    #will convert sparse (and dense matrixes) to the csc
    #format required by MDV
    if not isinstance(matrix, scipy.sparse.csc_matrix):
        matrix = scipy.sparse.csc_matrix(matrix)
   
    #check indexes are in sync and if so just return the csc matrix
    if list(main_names) == list(mod_names):
        return matrix,True
    #create lookup of mod indices to main indices
    main_map  =  {name: i for i, name in enumerate(main_names)}
    lookup = np.array([main_map[name] for name in mod_names])
    # Apply the lookup to the entire mod indices array using vectorized approach
    indices  = lookup[matrix.indices]
    # create a new sparse matrix 
    matrix = scipy.sparse.csc_matrix((matrix.data, indices, matrix.indptr), 
                                          shape=(len(main_names), matrix.shape[1]),
                                          dtype=matrix.dtype)
    return matrix,True



def convert_vcf_to_df(vcf_filename: str) -> pd.DataFrame:
    f = open(vcf_filename, "r")
    metadata = {}
    while True:
        line = f.readline()
        if line.startswith("##"):
            key, value = line[2:].strip().split("=", 1)
            if key in metadata:
                metadata[key].append(value)
            else:
                metadata[key] = [value]
        if line.startswith("#CHROM"):
            break
    # todo - do something with metadata
    # print(metadata)
    temp_file = "temp.csv"
    ## todo with tempfile
    with open(temp_file, "w") as tmp:
        while line:
            tmp.write(line)
            line = f.readline()
    f.close()

    df = pd.read_csv(temp_file, sep="\t")
    os.remove(temp_file)
    return df


def compute_vcf_end(df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute the end position of the variant determined from 'POS', 'REF' and 'ALT'.

    This is added as a column 'END' in the given DataFrame.
    """

    def varlen(s) -> int:
        return max([len(v) for v in str(s).split(",")])

    df["END"] = df["POS"] + df[["REF", "ALT"]].map(varlen).max(axis=1)
    return df


def convert_vcf_to_mdv(folder: str, vcf_filename: str) -> MDVProject:
    """
    Converts a VCF file to an MDV project.
    The VCF file must be tab-delimited, with the header lines starting with "##" and
    column names in the line starting with "#CHROM".

    An 'END' column is derived, which is the end position of the variant determined from 'POS', 'REF' and 'ALT'.
    """
    p = MDVProject(folder, delete_existing=True)
    df = convert_vcf_to_df(vcf_filename)
    # for single nucleotide variants, we still need an end position for the genome browser
    # maybe add_genome_browser should be able to understand only one position?
    # other forms VCF variants have a length, so we could use that...
    df = compute_vcf_end(df)
    # ^^ I should verify that this makes sense from biological perspective
    columns = [
        {"name": "#CHROM", "datatype": "text"},  # chromosome
        {"name": "POS", "datatype": "integer"},  # start of the variant
        {
            "name": "END",
            "datatype": "integer",
        },  # not standard VCF, but useful for genome browser(? maybe temporary)
        {
            "name": "ID",
            "datatype": "unique",
            "delimiter": ";",
        },  # should be unique, but also semicolon-delimited - could be useful to preserve this
        {"name": "REF", "datatype": "multitext", "delimiter": ","},  # reference base(s)
        {
            "name": "ALT",
            "datatype": "multitext",
            "delimiter": ",",
        },  # comma-delimited list of alternate non-reference alleles
        {
            "name": "QUAL",
            "datatype": "double",
        },  # phred-scaled quality score for the assertion made in ALT
        {
            "name": "FILTER",
            "datatype": "multitext",
            "delimiter": ";",
        },  # PASS or semicolon-delimited list of filters that the variant has failed
        {
            "name": "INFO",
            "datatype": "multitext",
            "delimiter": ",",
        },  # comma-delimited list of additional information, no white space, semicolons or equals signs permitted
        # ^^^ note, the first random vcf file I found has all manner of = and ; in the INFO field, so I'm not enclined to parse this too rigidly
        # {'name': 'FORMAT', 'datatype': 'text'}, # not sure why copilot thought this should be here - not in the VCF spec
    ]
    name = os.path.basename(vcf_filename)
    p.add_datasource(name, df, columns=columns)
    p.add_genome_browser(name, ["#CHROM", "POS", "END"])
    p.set_editable(True)
    return p


def create_regulamentary_project_from_pipeline(
        output,
        config,
        results_folder,
        atac_bw=None,
        peaks="merge",
        genome="hg38",
        openchrom="DNase"
):
    """Creates a regulamentary project from pipeline outputs.

    Args:
        output (str): Path to the directory which will house the MDV Project
        config (str): Path to the YAML configuration file.
        results_folder (str): Base path to the results directory.
        atac_bw (str, optional): Path to ATAC-seq bigWig file. Defaults to None.
        peaks (str, optional): Name of the peaks subdirectory. Defaults to "merge".
        genome (str, optional): Genome assembly version to use. Defaults to "hg38".
        openchrom (str, optional): Name of the open chromatin mark. Defaults to "DNase".

    Returns:
        An MDVProject 
    """
    fold = join(results_folder, peaks)
    marks = ["H3K4me1", "H3K4me3", "H3K27ac", "CTCF", "ATAC"]

    # Load configuration YAML
    with open(config, 'r') as file:
        info = yaml.safe_load(file)
    # was this created from bigwigs alone?
    from_bw = True if info.get("bigwigs") else False
    # Get the bed files for each mark
    # beds have been created in the pipeline
    if from_bw and info["bigwigs"].get("create_bed_files"):
        beds= {mark:join(results_folder,"bed_files",f"{mark}.bed") for mark in marks}
    #beds have been specified in the config
    else:
        beds = {mark: info["union_peaks"].get(f"bed_{mark}") for mark in marks}

    # Get the bigWig files for each mark
    if from_bw:
        #bigwigs including ATAC bigwig  specified in 'bigwigs' entry
        bigwigs={mark: info["bigwigs"].get(mark) for mark in marks}
    else:
        #otherwise bigwigs are specified in the 'compute_matrix_bigwigs' entry
        bigwigs = {mark: info["compute_matrix_bigwigs"].get(f"bigwig_{mark}") for mark in marks}
        bigwigs["ATAC"] = atac_bw  # Override ATAC with provided file if given

    # Set the path to the regulatory table
    table = join(fold, "08_REgulamentary", "mlv_REgulamentary.csv")

    # Determine genome version, considering blacklist genome override
    gen = genome
    bl = info.get("remove_blacklist")
    if bl and bl.get("genome"):
        gen = bl.get("genome")

    # Define matrix data source and region order
    matrix = {
        "data": join(fold, "09_metaplot", "matrix.csv"),
        "order": join(fold, "04_sort_regions", "sort_union.bed"),
        "marks": ["H3K4me1", "H3K4me3", "H3K27ac", "CTCF"]
    }

    return create_regulamentary_project(
        output,
        table,
        bigwigs,
        beds,
        matrix,
        openchrom=openchrom,
        genome=gen
    )

def create_regulamentary_project(
    output: str,
    table,
    bigwigs,
    beds,
    matrix=None,
    openchrom="DNase",
    marks = None,
    mark_colors = None,
    genome="hg38"
):
    """
    Creates a regulatory project visualization from input data sources.

    This method constructs a project using signal and peak files for 
    various histone marks and chromatin accessibility, adds them as data 
    sources and tracks, and configures a genome browser and visualization views.

    Args:
        output (str): Output directory or file for the project.
        table (str): Path to the CSV table containing regulatory element data.
        bigwigs (dict): Dictionary mapping mark names to bigWig file paths or URLs.
        beds (dict): Dictionary mapping mark names to BED file paths.
        matrix (dict or None, optional): Matrix and order file information for heatmaps, or None.
        openchrom (str, optional): Name for open chromatin mark. Defaults to "DNase".
        marks (list of str, optional): List of marks to process. Defaults to `["ATAC", "H3K4me1", "H3K4me3", "H3K27ac", "CTCF"]`.
        mark_colors (list of str, optional): List of colors for the marks. Defaults to a preset palette.
        genome (str, optional): Genome assembly to use. Defaults to "hg38".

    Returns:
        MDVProject: The project object constructed with the given data and views.

    """
    if marks is None:
        marks = ["ATAC", "H3K4me1", "H3K4me3", "H3K27ac", "CTCF"]
    if mark_colors is None:
        mark_colors = ["#eb9234", "#349beb", "#3aeb34", "#c4c41f", "#ab321a"]
    # Get the template directory
    tdir = join(split(os.path.abspath(__file__))[0], "templates")
    p = MDVProject(output, delete_existing=True)

    # Load regulatory elements table
    mdv = pd.read_csv(table, sep="\t")
    columns = [{"name": "start", "datatype": "int32"}, {"name": "end", "datatype": "int32"}]
    p.add_datasource("elements", mdv, columns)

    default_tracks = []
    for mark, color in zip(marks, mark_colors):
        name = mark if mark != 'ATAC' else openchrom
        bw = bigwigs.get(mark)
        if bw:
            url = bw
            if not bw.startswith("http"):
                fname = split(bw)[1]
                shutil.copy(bw, join(p.trackfolder, fname))
                url = f"./tracks/{fname}"
            default_tracks.append(
                {
                    "url": url,
                    "short_label": f"{name} cov",
                    "color": color,
                    "height": 60,
                    "track_id": f"coverage_{mark}"
                }
            )
        bed = beds.get(mark)
        if bed:
            url = bed
            if not bed.startswith("http"):
                fname = split(bed)[1]
                # Process BED files for browser compatibility
                if bed.endswith(".bed"):
                    # Bed file processing: remove header and keep first 3 columns
                    df = pd.read_csv(bed, sep="\t", header=None)
                    first_row = df.iloc[0]
                    cell_str  = str(first_row[1]).strip()
                    has_header = not cell_str.isdigit()
                    if has_header:
                        df = df.iloc[1:]
                    df = df.iloc[:, :3]
                    t_file = join(p.trackfolder, f"{fname}.temp")
                    o_file = join(p.trackfolder, fname)
                    df.to_csv(t_file, sep="\t", header=False, index=False)
                    create_bed_gz_file(t_file, o_file)
                    os.remove(t_file)
                    url = f"./tracks/{fname}.gz"
                else:
                    to_file = join(p.trackfolder, fname)
                    shutil.copyfile(bed, to_file)
                    # Copy tabix index if present
                    if bed.endswith(".gz"):
                        shutil.copyfile(f"{bed}.tbi", f"{to_file}.tbi")
                    url = f"./tracks/{fname}"
            default_tracks.append(
                {
                    "url": url,
                    "short_label": f"{name} peaks",
                    "color": color,
                    "height": 15,
                    "featureHeight": 5,
                    "track_id": f"peaks_{mark}"
                }
            )

    extra_params = {
        "default_tracks": default_tracks,
        "default_parameters": {
            "view_margins": {"type": "fixed_length", "value": 4000},
            "color_by": "RE",
            "color_legend": {"display": False, "pos": [5, 5]},
            "feature_label": "RE",
        },
    }
    p.add_genome_browser(
        "elements", ["chromosome", "start", "end"], extra_params=extra_params
    )
    p.add_refseq_track("elements", genome)
    if matrix:
        _create_dt_heatmap(p, matrix)

    # Load and extend visualization views
    with open(join(tdir, "views", "regulamentary.json")) as f:
        view = json.load(f)
    if not matrix:
        view = [x for x in view if x["type"] != "deeptools_heatmap"]
    gb = p.get_genome_browser("elements")
    gb.update(
        {
            "id": "browser",
            "size": [389, 550],
            "position": [1046, 8],
            "title": "elements",
        }
    )
    view.append(gb)
    p.set_view(
        "default",
        {
            "initialCharts": {
                "elements": view,
            }
        },
    )
    return p

def _create_dt_heatmap(
        project,
        matrix,
        ds="elements",
        ifields = None
):
    if ifields is None:
        ifields = ["chromosome", "start", "end"]
    # get the order of the regions in the matrix
    data = {x:project.get_column(ds,x) for x in ifields}
    mdv_i = pd.DataFrame(data)
    mdv_i = mdv_i.set_index(ifields)
    order = pd.read_csv(matrix["order"], sep="\t", header=None)
    order = order.set_index([0, 1, 2])
    order = order[order.index.isin(mdv_i.index)]
    mdv_i = mdv_i.reset_index()
    mdv_i["row_position"] = mdv_i.index
    arr = order.index.map(
        mdv_i.set_index(["chromosome", "start", "end"])["row_position"]
    )
    arr = np.array(arr, dtype=np.int32)

    # add the binary data
    output = join(project.dir, "binarydata", "elements")
    os.makedirs(output, exist_ok=True)
    hm = gzip.open(join(output, "heat_map.gz"), "wb")
    hm.write(arr.tobytes())

    # get the matrix and drop last two columns
    mt = pd.read_csv(matrix["data"], sep="\t")
    mt = mt.iloc[:, :1600]
    # flatten and halve in size (mean of adjacent values)
    arr = mt.values.flatten()
    arr = arr.reshape(-1, 2)
    arr = np.mean(arr, axis=1)
    # normalize and clip to 255 (1 byte per value)
    mx = np.percentile(arr, 99.99)
    arr = (arr / mx) * 255
    arr = np.clip(arr, 0, 255)
    arr = arr.astype(np.uint8)
    hm.write(arr.tobytes())
    hm.close()
    # add the metadata
    md = project.get_datasource_metadata(ds)
    md["deeptools"] = {
        "maps": {
            "default": {
                "data": "heat_map",
                "rows": md["size"],
                "cols": 800,
                "groups": matrix["marks"],
                "max_scale": mx,
            }
        }
    }
    md["binary_data_loader"] = True
    project.set_datasource_metadata(md)


def _add_dims(table, dims, max_dims: int, stub: str = "") -> pd.DataFrame:
    """
    Add dimensions to a table from a dictionary of dimension reductions.

    Args:
        table (pd.DataFrame): The table to add dimensions to.
        dims (dict): A dictionary of dimension reductions.
        max_dims (int): The maximum number of dimensions to add.
        stub (str): A string to prefix the dimension names with.

    Returns:
        pd.DataFrame: The table with dimensions added.
    """
    if len(dims.keys()) == 0:
        return table
    #stub is if there is  more than one modality e.g.
    #RNA UMAP, ATAC UMAP
    if stub:
        stub=stub+":"
    for dname,data in dims.items():
        #not sure whats going on here but cannot extract dims
        if len(data.shape) == 1:
            continue
        num_dims= min(max_dims,data.shape[1])
        #name the columns
        cols  = [f"{stub}{dname}_{i+1}" for i in range(num_dims)]
        dim_table= pd.DataFrame()
        #can either be an ndarray or a dataframe
        if isinstance(data, np.ndarray):
            dim_table = pd.DataFrame(data[:, 0:num_dims])
            dim_table.index = dims.dim_names
        elif isinstance(data, pd.DataFrame):
            # Should already have correct index
            dim_table = data.iloc[:, 0:num_dims]
        else:
            print (f"unrecognized dimension reduction format - {type(data)} for dim reduction {dname} ")
            continue
        dim_table.columns=cols
        #merge the tables to ensure the dims are in sync with the main table
        #perhaps only necessary with mudata objects but will do no harm
        table = table.merge(dim_table, left_index=True, right_index=True, how= "left")
    return table


def _default_merge_prefix(extra_dir: str) -> str:
    """
    Generate a deterministic prefix derived from the extra project directory.
    Ensures collisions are unlikely when importing datasources.
    """
    candidate = os.path.basename(os.path.normpath(extra_dir)) or "merged"
    candidate = re.sub(r"[^0-9A-Za-z_]+", "_", candidate)
    candidate = candidate.strip("_") or "merged"
    if not candidate.endswith("_"):
        candidate += "_"
    return candidate


def _replace_path_segment(path: str, segment: str, replacement: str) -> str:
    """
    Replace whole path segments inside a relative path.
    """
    if not isinstance(path, str):
        return path
    parts = path.split("/")
    changed = False
    for idx, part in enumerate(parts):
        if part == segment:
            parts[idx] = replacement
            changed = True
    return "/".join(parts) if changed else path


def merge_projects(
    base_dir: str,
    extra_dir: str,
    prefix: str | None = None,
    view_prefix: str | None = None,
) -> MDVProject:
    """
    Merge an existing MDV project into another project.

    Args:
        base_dir: Path to the target MDV project that will receive new datasources.
        extra_dir: Path to the MDV project to merge into ``base_dir``.
        prefix: Optional prefix applied to imported datasource names. When omitted a
            prefix derived from ``extra_dir`` is used.
        view_prefix: Optional prefix for imported views. Defaults to ``prefix`` with
            the trailing underscore trimmed.

    Returns:
        The updated ``MDVProject`` instance for ``base_dir``.
    """
    base_dir = os.path.abspath(base_dir)
    extra_dir = os.path.abspath(extra_dir)
    if base_dir == extra_dir:
        raise ValueError("base_dir and extra_dir must refer to different projects.")

    required_files = ("datafile.h5", "datasources.json", "views.json", "state.json")
    for directory in (base_dir, extra_dir):
        for fname in required_files:
            if not os.path.exists(os.path.join(directory, fname)):
                raise FileNotFoundError(
                    f"'{directory}' does not appear to be an MDV project (missing {fname})."
                )

    base_project = MDVProject(base_dir)
    extra_project = MDVProject(extra_dir)

    prefix_value = prefix
    if prefix_value is None:
        prefix_value = _default_merge_prefix(extra_dir)
    # Allow explicitly disabling prefixing via an empty string.
    prefix_value = prefix_value or ""

    base_names = set(base_project.get_datasource_names())
    extra_datasources = copy.deepcopy(extra_project.datasources)

    name_mapping: Dict[str, str] = {}
    for ds in extra_datasources:
        old_name = ds["name"]
        new_name = f"{prefix_value}{old_name}" if prefix_value else old_name
        if new_name in base_names:
            raise ValueError(
                f"Datasource '{new_name}' already exists in base project. "
                f"Provide a prefix to avoid collisions."
            )
        name_mapping[old_name] = new_name
        base_names.add(new_name)

    view_prefix_value = view_prefix
    if view_prefix_value is None:
        if prefix_value:
            view_prefix_value = prefix_value.rstrip("_")
        else:
            view_prefix_value = _default_merge_prefix(extra_dir).rstrip("_")
    if view_prefix_value:
        if not view_prefix_value.endswith("_"):
            view_prefix_value = f"{view_prefix_value}_"

    dir_tasks: List[Tuple[str, str]] = []
    file_tasks: List[Tuple[str, str]] = []
    track_tasks: List[Tuple[str, str]] = []
    region_field_updates: List[Tuple[str, str, Dict[str, str]]] = []
    string_replacements: Dict[str, str] = {k: v for k, v in name_mapping.items()}

    dir_destinations: set[str] = set()
    file_destinations: set[str] = set()

    def schedule_directory(src: str, dest: str) -> None:
        if not os.path.isdir(src):
            return
        if dest in dir_destinations:
            return
        dir_destinations.add(dest)
        dir_tasks.append((src, dest))

    def schedule_file(src: str, dest: str) -> None:
        if not os.path.isfile(src):
            return
        if dest in file_destinations:
            return
        file_destinations.add(dest)
        file_tasks.append((src, dest))

    def schedule_track(old_rel: str, new_rel: str) -> None:
        track_tasks.append((old_rel, new_rel))

    prepared_datasources: List[dict] = []

    def rewrite_json_path(value: str) -> str:
        if not isinstance(value, str) or not value.startswith("json/"):
            return value
        filename = value.split("/", 1)[1]
        new_filename = f"{prefix_value}{filename}" if prefix_value else filename
        new_rel = f"json/{new_filename}"
        schedule_file(os.path.join(extra_dir, value), os.path.join(base_dir, new_rel))
        if new_rel != value:
            string_replacements[value] = new_rel
        return new_rel

    for ds in extra_datasources:
        old_name = ds["name"]
        new_name = name_mapping[old_name]
        ds_meta = copy.deepcopy(ds)
        ds_meta["name"] = new_name

        links = ds_meta.get("links")
        if links:
            new_links = {}
            for target, link in links.items():
                new_links[name_mapping.get(target, target)] = link
            ds_meta["links"] = new_links

        # Binary data and datasource images
        schedule_directory(
            os.path.join(extra_dir, "binarydata", old_name),
            os.path.join(base_dir, "binarydata", new_name),
        )
        schedule_directory(
            os.path.join(extra_dir, "images", old_name),
            os.path.join(base_dir, "images", new_name),
        )

        images_meta = ds_meta.get("images")
        if isinstance(images_meta, dict):
            for img_name, img_info in images_meta.items():
                if not isinstance(img_info, dict):
                    continue
                base_url = img_info.get("base_url")
                if isinstance(base_url, str):
                    updated_url = _replace_path_segment(base_url, old_name, new_name)
                    if updated_url != base_url:
                        string_replacements[base_url] = updated_url
                        img_info["base_url"] = updated_url

        # Genome browser assets
        genome_browser = ds_meta.get("genome_browser")
        if isinstance(genome_browser, dict):
            sanitized_old = secure_filename(old_name)
            sanitized_new = secure_filename(new_name)

            def adjust_track(url: str | None) -> str | None:
                if not isinstance(url, str) or not url.startswith("tracks/"):
                    return url
                old_rel = url[len("tracks/") :]
                new_rel = old_rel
                if sanitized_old and sanitized_old != sanitized_new and old_rel.startswith(
                    sanitized_old
                ):
                    new_rel = sanitized_new + old_rel[len(sanitized_old) :]
                new_url = f"tracks/{new_rel}"
                schedule_track(old_rel, new_rel)
                if new_url != url:
                    string_replacements[url] = new_url
                return new_url

            default_track = genome_browser.get("default_track")
            if isinstance(default_track, dict):
                default_track["url"] = adjust_track(default_track.get("url"))
                if default_track.get("label") == old_name:
                    default_track["label"] = new_name

            default_tracks = genome_browser.get("default_tracks")
            if isinstance(default_tracks, list):
                for track in default_tracks:
                    if not isinstance(track, dict):
                        continue
                    track["url"] = adjust_track(track.get("url"))
                    if track.get("short_label") == old_name:
                        track["short_label"] = new_name

            atac_track = genome_browser.get("atac_bam_track")
            if isinstance(atac_track, dict):
                atac_track["url"] = adjust_track(atac_track.get("url"))

        # Regions metadata (spatial projects)
        regions = ds_meta.get("regions")
        if isinstance(regions, dict):
            all_regions = regions.get("all_regions", {})
            if isinstance(all_regions, dict) and all_regions:
                new_all_regions: Dict[str, dict] = {}
                region_mapping: Dict[str, str] = {}
                for region_name, region_info in all_regions.items():
                    if not isinstance(region_info, dict):
                        new_all_regions[region_name] = region_info
                        continue
                    new_region_name = (
                        f"{prefix_value}{region_name}" if prefix_value else region_name
                    )
                    updated_info = copy.deepcopy(region_info)
                    if (
                        isinstance(updated_info.get("spatial"), dict)
                        and updated_info["spatial"].get("file") == region_name
                    ):
                        updated_info["spatial"]["file"] = new_region_name
                    viv_image = updated_info.get("viv_image")
                    if isinstance(viv_image, dict):
                        for key in ("file", "url"):
                            if key in viv_image and isinstance(viv_image[key], str):
                                new_val = _replace_path_segment(
                                    viv_image[key], region_name, new_region_name
                                )
                                if new_val != viv_image[key]:
                                    string_replacements[viv_image[key]] = new_val
                                    viv_image[key] = new_val
                    images_dict = updated_info.get("images")
                    if isinstance(images_dict, dict):
                        for image_meta in images_dict.values():
                            if not isinstance(image_meta, dict):
                                continue
                            for key in ("url", "file"):
                                if key in image_meta and isinstance(image_meta[key], str):
                                    new_val = _replace_path_segment(
                                        image_meta[key], region_name, new_region_name
                                    )
                                    if new_val != image_meta[key]:
                                        string_replacements[image_meta[key]] = new_val
                                        image_meta[key] = new_val
                    json_path = updated_info.get("json")
                    if isinstance(json_path, str):
                        updated_info["json"] = rewrite_json_path(json_path)
                    new_all_regions[new_region_name] = updated_info
                    schedule_directory(
                        os.path.join(extra_dir, "spatial", region_name),
                        os.path.join(base_dir, "spatial", new_region_name),
                    )
                    if new_region_name != region_name:
                        region_mapping[region_name] = new_region_name
                        string_replacements[region_name] = new_region_name
                        string_replacements[f"spatial/{region_name}"] = (
                            f"spatial/{new_region_name}"
                        )
                regions["all_regions"] = new_all_regions
                if region_mapping:
                    region_field = regions.get("region_field")
                    if region_field:
                        region_field_updates.append(
                            (new_name, region_field, region_mapping)
                        )

        prepared_datasources.append(ds_meta)

    # Copy HDF5 groups
    with h5py.File(extra_project.h5file, "r") as extra_h5, h5py.File(
        base_project.h5file, "a"
    ) as base_h5:
        for old_name, new_name in name_mapping.items():
            if old_name not in extra_h5:
                raise KeyError(
                    f"Datasource '{old_name}' missing from extra project datafile."
                )
            if new_name in base_h5:
                raise ValueError(
                    f"HDF5 group '{new_name}' already exists in base project datafile."
                )
            extra_h5.copy(extra_h5[old_name], base_h5, name=new_name)

    # Update region field values within copied datasources (if required)
    if region_field_updates:
        with h5py.File(base_project.h5file, "a") as base_h5:
            for ds_name, field, region_map in region_field_updates:
                if not region_map:
                    continue
                group = base_h5.get(ds_name)
                if group is None or field not in group:
                    continue
                dataset = group[field]
                data = dataset[:]
                if data.dtype.kind == "S":
                    updated = np.array(
                        [
                            region_map.get(val.decode("utf-8"), val.decode("utf-8")).encode(
                                "utf-8"
                            )
                            for val in data
                        ],
                        dtype=data.dtype,
                    )
                    dataset[...] = updated
                elif data.dtype.kind == "U":
                    updated = np.array(
                        [region_map.get(val, val) for val in data], dtype=data.dtype
                    )
                    dataset[...] = updated
                elif data.dtype.kind == "O":
                    updated_values = []
                    for val in data:
                        if isinstance(val, bytes):
                            decoded = val.decode("utf-8")
                            updated_values.append(
                                region_map.get(decoded, decoded).encode("utf-8")
                            )
                        else:
                            decoded = str(val)
                            replacement = region_map.get(decoded, val)
                            updated_values.append(replacement)
                    dataset[...] = updated_values

    # Persist datasource metadata now that the HDF5 groups exist
    for ds_meta in prepared_datasources:
        base_project.set_datasource_metadata(ds_meta)

    # Convert track tasks into file copies (including .tbi siblings)
    for old_rel, new_rel in track_tasks:
        src = os.path.join(extra_dir, "tracks", old_rel)
        dest = os.path.join(base_dir, "tracks", new_rel)
        schedule_file(src, dest)
        src_tbi = f"{src}.tbi"
        dest_tbi = f"{dest}.tbi"
        if os.path.exists(src_tbi):
            schedule_file(src_tbi, dest_tbi)

    # Ensure target directories exist before copying.
    for subdir in ("binarydata", "images", "tracks", "spatial", "json"):
        Path(os.path.join(base_dir, subdir)).mkdir(parents=True, exist_ok=True)

    for src, dest in dir_tasks:
        if not os.path.isdir(src):
            continue
        if os.path.exists(dest):
            raise FileExistsError(
                f"Cannot merge '{src}' into '{dest}': destination already exists."
            )
        Path(dest).parent.mkdir(parents=True, exist_ok=True)
        shutil.copytree(src, dest)

    for src, dest in file_tasks:
        if not os.path.isfile(src):
            continue
        if os.path.exists(dest):
            raise FileExistsError(
                f"Cannot merge '{src}' into '{dest}': destination already exists."
            )
        Path(dest).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dest)

    # Merge views
    extra_views = copy.deepcopy(extra_project.views)

    def apply_replacements(obj):
        if isinstance(obj, dict):
            return {key: apply_replacements(value) for key, value in obj.items()}
        if isinstance(obj, list):
            return [apply_replacements(item) for item in obj]
        if isinstance(obj, str):
            return string_replacements.get(obj, obj)
        return obj

    def remap_view(view_config: dict) -> dict:
        view_copy = copy.deepcopy(view_config)
        if "initialCharts" in view_copy and isinstance(view_copy["initialCharts"], dict):
            remapped_initial = {}
            for ds_name, charts in view_copy["initialCharts"].items():
                new_key = string_replacements.get(ds_name, ds_name)
                remapped_initial[new_key] = [
                    apply_replacements(chart) for chart in charts
                ]
            view_copy["initialCharts"] = remapped_initial
        if "dataSources" in view_copy and isinstance(view_copy["dataSources"], dict):
            remapped_sources = {}
            for ds_name, config in view_copy["dataSources"].items():
                new_key = string_replacements.get(ds_name, ds_name)
                remapped_sources[new_key] = apply_replacements(config)
            view_copy["dataSources"] = remapped_sources
        return apply_replacements(view_copy)

    for view_name, view_data in extra_views.items():
        new_view_name = view_name
        if view_prefix_value:
            new_view_name = f"{view_prefix_value}{view_name}"
        if base_project.get_view(new_view_name):
            raise ValueError(
                f"View '{new_view_name}' already exists in base project. "
                f"Provide a different view prefix."
            )
        base_project.set_view(new_view_name, remap_view(view_data))

    logger.info(
        "Merged %d datasources from '%s' into '%s'.",
        len(prepared_datasources),
        extra_dir,
        base_dir,
    )
    return base_project
