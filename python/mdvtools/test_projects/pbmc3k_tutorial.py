## Before running this script:
# Run the following code to create the necessary folders and download the data
# Firstly cd into the folder containing the script
# Then uncomment and run the following from your terminal
# !mkdir data
# !wget http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz -O data/pbmc3k_filtered_gene_bc_matrices.tar.gz
# !cd data; tar -xzf pbmc3k_filtered_gene_bc_matrices.tar.gz
# !mkdir write



# https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html

import numpy as np
import pandas as pd
import scanpy as sc 
import json
from mdvtools.mdvproject import MDVProject
import os
# Import the chart classes
from mdvtools.charts.row_chart import RowChart
from mdvtools.charts.dot_plot import DotPlot
from mdvtools.charts.histogram_plot import HistogramPlot
from mdvtools.charts.scatter_plot_3D import ScatterPlot3D
from mdvtools.charts.scatter_plot import ScatterPlot
from mdvtools.charts.table_plot import TablePlot


# * * * Data Analysis section * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()

dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(dir)

adata = sc.read_10x_mtx(
    "data/filtered_gene_bc_matrices/hg19/",  # the directory with the `.mtx` file
    var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
    cache=True,  # write a cache file for faster subsequent reading
)

# adata.var_names_make_unique()  # this is unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx`

# annotate the group of mitochondrial genes as "mt"
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
)

# pylance in vscode is giving type errors for the following lines - even though the code not only runs correctly,
# but pyright when called from terminal or in CI action also doesn't give any errors.
adata = adata[adata.obs.n_genes_by_counts < 2500, :] # type: ignore
adata = adata[adata.obs.pct_counts_mt < 5, :].copy() # type: ignore

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

adata.raw = adata

adata = adata[:, adata.var.highly_variable] # type: ignore

sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])

sc.pp.scale(adata, max_value=10)

sc.tl.pca(adata, svd_solver="arpack")

sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)

sc.tl.leiden(
    adata,
    resolution=0.9,
    random_state=0,
    #flavor="igraph",
    n_iterations=2,
    directed=False,
)

sc.tl.rank_genes_groups(adata, "leiden", method="t-test")
sc.settings.verbosity = 2  # reduce the verbosity
sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
sc.tl.rank_genes_groups(adata, "leiden", method="logreg", max_iter=1000)

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

cells_df = pd.DataFrame(adata.obs)

pca_np = np.array(adata.obsm["X_pca"])
cells_df["X_pca_1"] = pca_np[:, 0]
cells_df["X_pca_2"] = pca_np[:, 1]
cells_df["X_pca_3"] = pca_np[:, 2]

umap_np = np.array(adata.obsm["X_umap"])
cells_df["X_umap_1"] = umap_np[:, 0]
cells_df["X_umap_2"] = umap_np[:, 1]

cells_df["cell_id"] = adata.obs.index
'''adata.layers['counts'] = adata.X.copy()

marker_genes = [
    *["IL7R", "CD79A", "MS4A1", "CD8A", "CD8B", "LYZ", "CD14"],
    *["LGALS3", "S100A8", "GNLY", "NKG7", "KLRB1"],
    *["FCGR3A", "MS4A7", "FCER1A", "CST3", "PPBP"],
]

#cells_df[marker_genes] = adata.var["gene_ids"][adata.var["gene_ids"] == marker_genes]

print(adata.uns["rank_genes_groups"]["names"].field(1))

result = adata.uns["rank_genes_groups"]
groups = result["names"].dtype.names
print(pd.DataFrame(
    {
        group + "_" + key[:1]: result[key][group]
        for group in groups
        for key in ["names", "scores"]
    }
)) '''
assert isinstance(adata.X, np.ndarray)
adata.layers['counts'] = adata.X.copy()

## these can be added dynamically with the 'Add Gene Expr' UI in the frontend
assert isinstance(adata.layers['counts'], np.ndarray)
counts = np.array(adata.layers['counts'])
transpose_counts = np.transpose(counts)
cells_df['ARVCF']=transpose_counts[adata.var.index.get_loc("ARVCF")]
#cells_df['CTB-113I20.2']=transpose_counts[adata.var.index.get_loc("CTB-113I20.2")])
cells_df['DOK3']=transpose_counts[adata.var.index.get_loc("DOK3")]
cells_df['FAM210B']=transpose_counts[adata.var.index.get_loc("FAM210B")]
cells_df['GBGT1']=transpose_counts[adata.var.index.get_loc("GBGT1")]
cells_df['NFE2L2']=transpose_counts[adata.var.index.get_loc("NFE2L2")]
#cells_df['PPBP']=transpose_counts[adata.var.index.get_loc("PPBP")])
cells_df['UBE2D4']=transpose_counts[adata.var.index.get_loc("UBE2D4")]
cells_df['YPEL2']=transpose_counts[adata.var.index.get_loc("YPEL2")]

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
#Αdd a second datasource 'genes'
gene_table = adata.var
gene_table["gene_ids"]=gene_table.index

varm_np = np.array(adata.varm["PCs"])
gene_table["PCs_1"] = varm_np[:, 0]
gene_table["PCs_2"] = varm_np[:, 1]

print(adata.var)
#print(gene_table.columns.values.tolist())

#gene_table["PCs_1"] = rna.varm["PCs"][:, 0]
#gene_table["PCs_2"] = rna.varm["PCs"][:, 1]

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
# CELLS GRAPHS
# Create a RowChart instance with configuration based on the JSON
row_chart = RowChart(
    title="leiden",
    param="leiden",
    position=[380, 270],
    size=[260, 260]
)

# Configure the row chart
row_chart.set_axis_properties("x", {"textSize": 13, "label": "", "tickfont": 10})

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

# Create a DotPlot instance with configuration based on the JSON
dot_plot = DotPlot(
    title="leiden",
    params=["leiden", "ARVCF", "DOK3", "FAM210B", "GBGT1", "NFE2L2", "UBE2D4", "YPEL2"],
    size=[400, 250],
    position=[10, 10]
)

# Configure the dot plot
dot_plot.set_axis_properties("x", {"label": "", "textSize": 13, "tickfont": 10})
dot_plot.set_axis_properties("y", {"label": "", "textSize": 13, "tickfont": 10})
dot_plot.set_axis_properties("ry", {"label": "", "textSize": 13, "tickfont": 10})
dot_plot.set_color_scale(log_scale=False)
dot_plot.set_color_legend(True, [40, 10])
dot_plot.set_fraction_legend(True, [0, 0])


# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
histogram_plot = HistogramPlot(
    title="n_genes_by_counts",
    param="n_genes_by_counts",
    bin_number=17,
    display_min=0,
    display_max=2500,
    size=[360, 250],
    position=[10, 280]
)

histogram_plot.set_x_axis(size=30, label="n_genes_by_counts", textsize=13, tickfont=10)
histogram_plot.set_y_axis(size=45, label="frequency", textsize=13, tickfont=10, rotate_labels=False)

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
# Create a ScatterPlot3D instance
scatter_plot = ScatterPlot3D(
    title="X_pca_1 x X_pca_2 x X_pca_3",
    params=["X_pca_1", "X_pca_2", "X_pca_3"],
    size=[250, 250],
    position=[420, 10],
    default_color="#377eb8",
    brush="default",
    center=[0, 0, 0],
    on_filter="hide",
    radius=5,
    opacity=0.8,
    axis_scales=[220, 220, 220],
    camera={"distance": 37543.999999999956, "theta": -1.038, "phi": 0.261}
    #color_by="leiden"
)

scatter_plot.set_color_by("leiden")

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
# Genes graphs
# Create a ScatterPlot instance
scatter2D_plot = ScatterPlot(
    title="PCs_1 x PCs_2",
    params=["PCs_1", "PCs_2"],
    size=[250, 250],
    position=[270, 10]
)
scatter2D_plot.set_default_color("#377eb8")
scatter2D_plot.set_brush("poly")
scatter2D_plot.set_opacity(0.8)
scatter2D_plot.set_radius(2)
scatter2D_plot.set_color_legend(display=True, position=[20, 1])
scatter2D_plot.set_axis_properties("x", {"label": "PCs_1", "textSize": 13, "tickfont": 10})
scatter2D_plot.set_axis_properties("y", {"label": "PCs_2", "textSize": 13, "tickfont": 10})
scatter2D_plot.set_color_by("dispersions")

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
table_plot = TablePlot(
    title="Data table",
    params= gene_table.columns.values.tolist() ,
    size=[250, 500],
    position=[10, 10]
)
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
'''histogram_plot = HistogramPlot(
    title="n_genes_by_counts",
    param="n_genes_by_counts",
    bin_number=17,
    display_min=0,
    display_max=2500,
    size=[360, 250],
    position=[10, 280]
)

histogram_plot.set_x_axis(size=30, label="n_genes_by_counts", textsize=13, tickfont=10)
histogram_plot.set_y_axis(size=45, label="frequency", textsize=13, tickfont=10, rotate_labels=False) '''

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

# Set up and serve the MDV project
base = os.path.expanduser('~/mdv')
project_path = os.path.join(base, 'pbmc3k')
p = MDVProject(os.path.expanduser(project_path), delete_existing=True)

# Add data to the project
p.add_datasource("cells", cells_df)
p.add_datasource("genes", gene_table)




# Convert plot data to JSON and setup the project view
list_charts_cells = []
list_charts_genes = []

row_chart_data = json.loads(json.dumps(row_chart.plot_data, indent=2).replace("\\", ""))
dot_chart_data = json.loads(json.dumps(dot_plot.plot_data, indent=2).replace("\\", ""))
histogram_chart_data = json.loads(json.dumps(histogram_plot.plot_data, indent=2).replace("\\", ""))
scatter_chart = json.loads(json.dumps(scatter_plot.plot_data, indent=2).replace("\\", ""))

table_chart = json.loads(json.dumps(table_plot.plot_data, indent=2).replace("\\", ""))
scatter2D_chart = json.loads(json.dumps(scatter2D_plot.plot_data, indent=2).replace("\\", ""))

list_charts_cells.append(row_chart_data)
list_charts_cells.append(dot_chart_data)
list_charts_cells.append(histogram_chart_data)
list_charts_cells.append(scatter_chart)

list_charts_genes.append(table_chart)
list_charts_genes.append(scatter2D_chart)

view_config = {'initialCharts': {"cells": list_charts_cells, "genes":list_charts_genes}}

#link the two datasets
p.add_rows_as_columns_link("cells","genes","gene_ids","Gene Expr")

#add the gene expression 
p.add_rows_as_columns_subgroup("cells","genes","gs",adata.X,name="scores",label="Gene Scores")

p.set_view("default", view_config)
p.set_editable(True)
p.serve()