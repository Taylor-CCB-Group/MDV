from scanpy import AnnData
import scipy
import pandas as pd
from os.path import join, split, exists
import os
import re
import json
import gzip
import copy
import yaml
import shutil
from math import ceil
import h5py
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Any
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
    #add drs derived from modalities (cell clustering)
    for name,mod in md.mod.items():
        table = _add_dims(table,mod.obsm,max_dims,name)
    table["cell_id"] = md.obs.index
    columns= [{"name":"cell_id","datatype":"unique"}]
    p.add_datasource("cells",table,columns)

    for mod in md.mod.keys():
        mdata = md.mod[mod]

        #adds the index to the data as a name column
        #This is usually the unique gene name - but not always
        if "name" in mdata.var.columns:
            logger.warning(f"name in modality {mod} will be overwritten with the index")
        #Add name column to the actual dataframe, previously it was added to the datasource
        #after creation. This way it ensures that the dataframe has a least one column, which
        #is required for subsequent steps
        mdata.var["name"]=mdata.var.index

        #add any DRs to the modality (e.g. gene/protein clustering)
        mdata.var = _add_dims(mdata.var,mdata.varm,max_dims)
     
        #add the modality to project
        p.add_datasource(mod,mdata.var)
        
        #mod is used as both the tag and the label
        #The name column is specified as the identifier that the user will use
        #It is derived from the index and is usually the gene 'name'
        #However it may not be appropriate and can be changed later on 
        p.add_rows_as_columns_link("cells",mod,"name",mod)
        matrix,sparse= get_matrix(mdata.X,md.obs_names,mdata.obs_names)
        #sometimes X is empty - all the data is in the layers
        if matrix.shape[1] !=0:
            p.add_rows_as_columns_subgroup("cells",mod,mod+"_expr",matrix,sparse=sparse, chunk_data=chunk_data)
        #now add the layers (if there are any)
        layers = mdata.layers.keys()
        for layer in layers:
            matrix  = mdata.layers[layer]
            matrix,sparse = get_matrix(matrix,md.obs_names,mdata.obs_names)
            p.add_rows_as_columns_subgroup("cells",mod,f"{mod}_{layer}",matrix,sparse=sparse, chunk_data=chunk_data)
    return p


# The main_names correspond to obs_names in the main anndata object
# and the mod_names those in the modality.
# Only required if data is being added from a modality,as sometimes
# the modality's obs_names will be in a different order and/or a subset of the main names
# Hence a sparse matrix corresponding to the main indices needs to be created
def get_matrix(matrix,main_names=None,mod_names=None) -> tuple[scipy.sparse.csc_matrix | np.ndarray, bool]:
    """
    Process and align a matrix to match main observation names.
    
    Args:
        matrix: Input matrix (sparse or dense, may be wrapped)
        main_names: Target observation names (e.g., from main AnnData object)
        mod_names: Source observation names (e.g., from modality)
    
    Returns:
        Tuple of (aligned_matrix, is_sparse)
    """
    if main_names is None:
        main_names = []
    if mod_names is None:
        mod_names = []
    
    # Unwrap matrix if needed
    matrix = matrix.value if hasattr(matrix,"value") else matrix
    
    # Handle backed sparse matrices
    if hasattr(matrix,"backend"):
        matrix = matrix._to_backed()
    
    # Validate matrix
    if matrix is None:
        raise ValueError("Matrix cannot be None")
    
    # Check for empty matrix
    if hasattr(matrix, 'shape') and (matrix.shape[0] == 0 or matrix.shape[1] == 0):
        logger.warning(f"Empty matrix with shape {matrix.shape}")
        return matrix, scipy.sparse.issparse(matrix)
    
    is_sparse = scipy.sparse.issparse(matrix)
    
    # Convert sparse matrix to CSC format for efficient column operations
    if is_sparse and not isinstance(matrix, scipy.sparse.csc_matrix):
        matrix = scipy.sparse.csc_matrix(matrix)
    
    # Early return if no alignment needed
    if  not any(main_names)  or not any(mod_names):
        return matrix, is_sparse
    
    # Validate dimensions
    if len(mod_names) != matrix.shape[0]:
        raise ValueError(
            f"Length mismatch: mod_names has {len(mod_names)} entries "
            f"but matrix has {matrix.shape[0]} rows"
        )
    
    # Quick check if already aligned
    if len(main_names) == len(mod_names) and all(m == n for m, n in zip(main_names, mod_names)):
        return matrix, is_sparse
    
    # Create lookup from main names to indices
    main_map = {name: i for i, name in enumerate(main_names)}
    
    # Identify missing and extra names
    missing_names = [name for name in mod_names if name not in main_map]
    
    if missing_names:
        n_missing = len(missing_names)
        if n_missing <= 5:
            logger.warning(
                f"{n_missing} cell(s) present in modality but absent in main object "
                f"and will be ignored: {', '.join(missing_names)}"
            )
        else:
            logger.warning(
                f"{n_missing} cells present in modality but absent in main object "
                f"and will be ignored (showing first 5): {', '.join(missing_names[:5])}"
            )
    
    # Check for duplicate names in mod_names
    if len(mod_names) != len(set(mod_names)):
        duplicates = [name for name in set(mod_names) if mod_names.count(name) > 1]
        logger.warning(
            f"Duplicate names found in modality: {', '.join(duplicates[:5])}. "
            f"Only the first occurrence will be used."
        )
    
    # Filter valid entries and build lookup
    valid_indices = []
    valid_mod_names = []
    seen = set()
    
    for i, name in enumerate(mod_names):
        if name in main_map and name not in seen:
            valid_indices.append(i)
            valid_mod_names.append(name)
            seen.add(name)
    
    if not valid_indices:
        raise ValueError(
            "No overlapping cell names between main object and modality. "
            "Cannot align matrices."
        )
    
    # Extract valid rows
    if len(valid_indices) < len(mod_names):
        valid_indices_array = np.array(valid_indices)
        matrix = matrix[valid_indices_array, :]
        mod_names = valid_mod_names
    
    # Build mapping from mod indices to main indices
    lookup = np.array([main_map[name] for name in mod_names], dtype=np.int32)
    
    if not is_sparse:
        # Dense matrix handling
        new_matrix = np.full((len(main_names), matrix.shape[1]), np.nan, dtype=np.float32)
        new_matrix[lookup, :] = matrix
        return new_matrix, False
    else:
      
        # Remap row indices to match main_names order
        new_indices = lookup[matrix.indices]
        
        # Verify indices are within bounds
        if np.any(new_indices < 0) or np.any(new_indices >= len(main_names)):
            raise ValueError("Index remapping produced out-of-bounds indices")
        
        # Create new sparse matrix with remapped rows
        try:
            new_matrix = scipy.sparse.csc_matrix(
                (matrix.data, new_indices, matrix.indptr),
                shape=(len(main_names), matrix.shape[1]),
                dtype=matrix.dtype
            )
        except Exception as e:
            logger.error(f"Failed to create remapped sparse matrix: {e}")
            raise
        
        return new_matrix, True



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

def create_capsequm_project(
    output: str,
    results_folder: str,
    extra_bigwigs: list =[]
):
    """
    Creates a CAPSeqUM project from the results folder.

    Args:
        output (str): Path to the output directory for the MDV project.
        results_folder (str): Base path to the results directory.
        extra_bigwigs (list, optional): List of additional bigWig files to include. Defaults to None.

    Returns:
        MDVProject: The created MDV project object.
    """
    p = MDVProject(output, delete_existing=True)
    # get info from config file
    with open(join(results_folder,"config.yml"), 'r') as file:
        config = yaml.safe_load(file)
    # read  in the oligo data
    df = pd.read_csv(join(results_folder, "all_oligos.tsv"), sep="\t")
    #add extra column merging is_best_oligo and pass_filter
    df["oligo_status"]=["best" if bool(x) else "pass" if bool(y) else "fail"\
                         for x,y in zip(df["is_best_oligo"], df["pass_filter"])]
    #get first best oligo to set genome browser location
    best_df = df[df["oligo_status"] == "best"]
    first_best = best_df.iloc[0] if not best_df.empty else df.iloc[0]
    chr_val = first_best["chr"]
    start_val = first_best["start"] -500
    stop_val = first_best["stop"] + 500

    #customize some of the columns
    cols = [
        {"name": "sequence", "datatype": "unique"},
        {
            "name": "oligo_status", 
            "datatype": "text",
            "values": ["best", "pass", "fail"],
            "colors": ["#4CAF50", "#9D9E48", "#F44336"],
        },
        {"name": "start", "datatype": "int32"},
        {"name": "stop", "datatype": "int32"}
    ]
    p.add_datasource("oligos", df, columns=cols)

    # make a datasource of the regions
    regions_file = join(results_folder, "regions", "regions.bed")
    #does it have features i.e. SNPs
    features_file  = join(results_folder, "regions","features.bed")
    has_features = exists(features_file)
  
    rds  = pd.read_csv(features_file if has_features else regions_file, sep="\t", header=None)
    rds.columns = ["chr", "start", "stop","region"]
    #get failed regions
    failed_file = join(results_folder, "failed_regions.txt")
    failed= set()
    if os.path.exists(failed_file):
        failed = pd.read_csv(failed_file, sep="\t", header=None)
        failed.columns = ["region"]
        failed = set(failed["region"].tolist()) 
    #add a column to the regions dataframe showing if oligos were found
    rds["oligos found"] = ["No" if x in failed else "Yes" for x in rds["region"]]
    p.add_datasource("regions", rds, columns=[ 
        {"name": "start", "datatype": "int32"},
        {"name": "stop", "datatype": "int32"},
        {"name": "region", "datatype": "unique"},
        {"name": "oligos found", "datatype": "text", 
         "values":["Yes","No"],"colors":["#4CAF50","#F44336"]}
    ])

    op = config["oligo_parameters"]
    # for features (snps) - need wide view margins to see all oligos
    vm = op["length"]*2 +50 if has_features else 500
    p.add_genome_browser(
        "regions", ["chr", "start", "stop"],name="all_regions",
        extra_params={
            "default_parameters": {
                "color_by": "oligos found",
                "highlight_selected_region": True,
                "view_margins": {"type": "fixed_length", "value": vm},
                "color_legend": {"display": False, "pos": [5, 5]}
            }
        }
    )

    #calculate height of track based on overlap i.e. how stacked the oligos are
    
    ot_height = ceil(op["length"]/op["step"]) * 7.5
    extra_params = {
        "default_parameters": {
            "color_by": "oligo_status",
            "highlight_selected_region":True,
            # will update when regions are selected
            "sync_with_datastores": [
                "regions"
            ],
            "view_margins": {"type": "fixed_length", "value": 500},
             "genome_location": {
                "chr": chr_val,
                "start": int(start_val),
                "end": int(stop_val)
            },
            "color_legend": {"display": False, "pos": [5, 5]}
        },
        "default_track_parameters":{
            "height":ot_height,
            "displayMode":"SQUISHED"
        }
    }

    p.add_genome_browser("oligos", ["chr", "start", "stop"],extra_params=extra_params)
    if config.get("genome"):
        try:
            p.add_refseq_track("oligos", config["genome"])
        except Exception as e:
           logger.warning(f"Could not add refseq track for genome '{config['genome']}'. Reason: {e}")

    tracks =[]
    if config.get("create_offtarget_bigwig"):
        tracks.append({"file": join(results_folder, "offtarget.bw"),"color":"#72211F",})
   
    
    if has_features:
        tracks.append({"file": features_file,"hideLabels":True,"color":"#6B81E4"})
    original = join(results_folder, "summits", "original_regions.bed")
    hideLabels=False
    if exists(original):
        tracks.append({"file": original,"color":"#A0A0A0" })
        hideLabels =True
    tracks.append({
        "file": join(results_folder, "regions", "regions.bed"),
        "hideLabels":hideLabels,
        "name":"Regions" if not exists(original) else "Summit Regions",
        "color": "#6B81E4",
    })
    for bws in extra_bigwigs:
        tracks.append({"file": bws,"color":"#463B86",})
    p.add_tracks("oligos", tracks )

    gb = p.get_genome_browser("oligos")
    gb.update(
        {
            "id": "browser",
            "size": [838,460],
            "position": [5,5],
            "title": "oligos",
        }
    )
    tdir = join(split(os.path.abspath(__file__))[0], "templates")
    with open(join(tdir, "views", "capsequm.json")) as f:
        views = json.load(f)
    
    views["oligos"].append(gb)
    p.set_view(
        "default",
        {
            "initialCharts": {
                "oligos": views["oligos"],
                "regions": views["regions"]
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
        if not isinstance(value, str):
            return value
        # Handle both json/ and images/ prefixes for geo.json files
        if value.startswith("json/"):
            prefix_dir = "json"
        elif value.startswith("images/"):
            prefix_dir = "images"
        else:
            return value
        filename = value.split("/", 1)[1]
        new_filename = f"{prefix_value}{filename}" if prefix_value else filename
        new_rel = f"{prefix_dir}/{new_filename}"
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
                # First pass: collect all actual spatial directory names from metadata
                # The actual spatial directories are named after sdata_name (e.g., "sample1.zarr"),
                # not the region names (e.g., "sample1_region1")
                spatial_dirs_to_copy: Dict[str, str] = {}  # old_spatial_dir -> new_spatial_dir
                
                for region_name, region_info in all_regions.items():
                    if not isinstance(region_info, dict):
                        continue
                    # Extract the actual spatial directory name from metadata
                    spatial_info = region_info.get("spatial")
                    if isinstance(spatial_info, dict):
                        spatial_file = spatial_info.get("file")
                        if isinstance(spatial_file, str):
                            # Track which spatial directories need to be copied
                            if spatial_file not in spatial_dirs_to_copy:
                                new_spatial_dir = (
                                    f"{prefix_value}{spatial_file}" if prefix_value else spatial_file
                                )
                                spatial_dirs_to_copy[spatial_file] = new_spatial_dir
                
                # Copy all spatial directories (not just those referenced by region names)
                # Also scan the actual spatial directory for any directories not in metadata
                extra_spatial_dir = os.path.join(extra_dir, "spatial")
                if os.path.isdir(extra_spatial_dir):
                    for item in os.listdir(extra_spatial_dir):
                        item_path = os.path.join(extra_spatial_dir, item)
                        if os.path.isdir(item_path):
                            # If not already scheduled, add it
                            if item not in spatial_dirs_to_copy:
                                new_spatial_dir = (
                                    f"{prefix_value}{item}" if prefix_value else item
                                )
                                spatial_dirs_to_copy[item] = new_spatial_dir
                
                # Schedule copying of all spatial directories
                for old_spatial_dir, new_spatial_dir in spatial_dirs_to_copy.items():
                    schedule_directory(
                        os.path.join(extra_dir, "spatial", old_spatial_dir),
                        os.path.join(base_dir, "spatial", new_spatial_dir),
                    )
                    if old_spatial_dir != new_spatial_dir:
                        string_replacements[f"spatial/{old_spatial_dir}"] = (
                            f"spatial/{new_spatial_dir}"
                        )
                        string_replacements[old_spatial_dir] = new_spatial_dir
                
                # Second pass: update region metadata with new names and paths
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
                    
                    # Update spatial file reference if it exists
                    spatial_info = updated_info.get("spatial")
                    if isinstance(spatial_info, dict):
                        spatial_file = spatial_info.get("file")
                        if isinstance(spatial_file, str):
                            # Update to new spatial directory name if it was renamed
                            new_spatial_file = spatial_dirs_to_copy.get(spatial_file, spatial_file)
                            if new_spatial_file != spatial_file:
                                spatial_info["file"] = new_spatial_file
                                string_replacements[spatial_file] = new_spatial_file
                    
                    # Update viv_image paths
                    viv_image = updated_info.get("viv_image")
                    if isinstance(viv_image, dict):
                        for key in ("file", "url"):
                            if key in viv_image and isinstance(viv_image[key], str):
                                # Update paths that reference old spatial directory names
                                old_val = viv_image[key]
                                new_val = old_val
                                for old_dir, new_dir in spatial_dirs_to_copy.items():
                                    if old_dir != new_dir and old_dir in old_val:
                                        new_val = old_val.replace(old_dir, new_dir)
                                        break
                                # Also update region name references
                                if new_val == old_val:
                                    new_val = _replace_path_segment(
                                        old_val, region_name, new_region_name
                                    )
                                if new_val != old_val:
                                    string_replacements[old_val] = new_val
                                    viv_image[key] = new_val
                    
                    # Update images dictionary paths
                    images_dict = updated_info.get("images")
                    if isinstance(images_dict, dict):
                        for image_meta in images_dict.values():
                            if not isinstance(image_meta, dict):
                                continue
                            for key in ("url", "file"):
                                if key in image_meta and isinstance(image_meta[key], str):
                                    old_val = image_meta[key]
                                    new_val = old_val
                                    # Update paths that reference old spatial directory names
                                    for old_dir, new_dir in spatial_dirs_to_copy.items():
                                        if old_dir != new_dir and old_dir in old_val:
                                            new_val = old_val.replace(old_dir, new_dir)
                                            break
                                    # Also update region name references
                                    if new_val == old_val:
                                        new_val = _replace_path_segment(
                                            old_val, region_name, new_region_name
                                        )
                                    if new_val != old_val:
                                        string_replacements[old_val] = new_val
                                        image_meta[key] = new_val
                    
                    json_path = updated_info.get("json")
                    if isinstance(json_path, str):
                        updated_info["json"] = rewrite_json_path(json_path)
                    
                    new_all_regions[new_region_name] = updated_info
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
                        # Also update the values array in the column metadata
                        columns = ds_meta.get("columns", [])
                        for col in columns:
                            if col.get("name") == region_field or col.get("field") == region_field:
                                if "values" in col and isinstance(col["values"], list):
                                    col["values"] = [
                                        region_mapping.get(v, v) for v in col["values"]
                                    ]
                                break

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
                if not isinstance(group, h5py.Group):
                    continue
                if field not in group:
                    continue
                dataset = group[field]
                if not isinstance(dataset, h5py.Dataset):
                    continue
                data = dataset[:]
                if data.dtype.kind == "S":
                    updated = np.array(
                        [
                            (
                                lambda v: (
                                    decoded := v.decode("utf-8"),
                                    mapped := region_map.get(decoded, decoded),
                                    mapped if isinstance(mapped, str) else decoded,
                                )[2]
                            )(val).encode("utf-8")
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
    
    # DEBUG
    print(f"DEBUG: extra_views keys: {list(extra_views.keys())}")
    print(f"DEBUG: view_prefix_value: '{view_prefix_value}'")

    def apply_replacements(obj: Any) -> Any:
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
        # DEBUG
        print(f"DEBUG: Added view '{new_view_name}'")

    logger.info(
        "Merged %d datasources from '%s' into '%s'.",
        len(prepared_datasources),
        extra_dir,
        base_dir,
    )
    return base_project
