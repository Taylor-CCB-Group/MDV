from scanpy import AnnData
import scipy
import pandas as pd
from os.path import join, split
import os
from .mdvproject import MDVProject,create_bed_gz_file
import numpy as np
import json
import gzip
import copy
import yaml
import shutil


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
def get_matrix(matrix,main_names=[],mod_names=[]):
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


def _add_dims(table, dims, max_dims,stub=""):
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
