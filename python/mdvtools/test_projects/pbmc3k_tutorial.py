import numpy as np
import pandas as pd
import scanpy as sc

import os
import json
from mdvtools.mdvproject import MDVProject
from mdvtools.charts.row_chart import RowChart
from mdvtools.charts.dot_plot import DotPlot
from mdvtools.charts.histogram_plot import HistogramPlot
from mdvtools.charts.scatter_plot_3D import ScatterPlot3D
from mdvtools.charts.scatter_plot import ScatterPlot
from mdvtools.charts.table_plot import TablePlot

# scanpy parameters for feedback level setting
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header() # printing a header of introductory information about the environment/library used

# reading the 3k 10x data

adata = sc.read_10x_mtx(
    "data/filtered_gene_bc_matrices/hg19/",  # the directory with the `.mtx` file
    var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
    cache=True,)  # write a cache file for faster subsequent reading)

# annotate the group of mitochondrial genes as "mt"
adata.var["mt"] = adata.var_names.str.startswith("MT-")

# running quality control metrics to filter the anndata object
sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

# using filter_cells and filter_genes to filter our cells
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# slicing the anndata object to remove genes that have too many mitochondrial genes expressed or too many total counts
adata = adata[adata.obs.n_genes_by_counts < 2500, :] 
adata = adata[adata.obs.pct_counts_mt < 5, :].copy()

# total-count normalising the anndata object to 10,000 reads per cell, so that counts become comparable among cells
sc.pp.normalize_total(adata, target_sum=1e4)

# logarithmising the data
sc.pp.log1p(adata)

# identifying highly-variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

# freezing the state of the anndata object to the normalised and logarithmised raw gene expression
adata.raw = adata

# filtering the highly-variable genes
adata = adata[:, adata.var.highly_variable]

# regressing out effects of total counts per cell and the percentage of mitochondrial genes expressed
sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])

# scaling the data to unit variance
sc.pp.scale(adata, max_value=10)

# running PCA
sc.tl.pca(adata, svd_solver="arpack")

# computing the neighborhood graph of cells using the PCA representation of the anndata object
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

sc.tl.umap(adata)

sc.tl.leiden(
    adata,
    resolution=0.9,
    random_state=0,
    #flavor="igraph", default is leidenalg
    n_iterations=2,
    directed=False,)

sc.tl.rank_genes_groups(adata, "leiden", method="t-test")
sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
sc.tl.rank_genes_groups(adata, "leiden", method="logreg", max_iter=1000)

# cells dataframe 

cells_df = pd.DataFrame(adata.obs)

# adding the PCA data to the cells dataframe
pca_np = np.array(adata.obsm["X_pca"])
cells_df["X_pca_1"] = pca_np[:, 0]
cells_df["X_pca_2"] = pca_np[:, 1]
cells_df["X_pca_3"] = pca_np[:, 2]

# adding the umap data to the cells dataframe
umap_np = np.array(adata.obsm["X_umap"])
cells_df["X_umap_1"] = umap_np[:, 0]
cells_df["X_umap_2"] = umap_np[:, 1]

cells_df["cell_id"] = adata.obs.index

assert isinstance(adata.X, np.ndarray)
adata.layers['counts'] = adata.X.copy()

## these can be added dynamically with the 'Add Gene Expr' UI in the frontend
assert isinstance(adata.layers['counts'], np.ndarray)
counts = np.array(adata.layers['counts'])
transpose_counts = np.transpose(counts)
cells_df['ARVCF']=transpose_counts[adata.var.index.get_loc("ARVCF")]
cells_df['DOK3']=transpose_counts[adata.var.index.get_loc("DOK3")]
cells_df['FAM210B']=transpose_counts[adata.var.index.get_loc("FAM210B")]
cells_df['GBGT1']=transpose_counts[adata.var.index.get_loc("GBGT1")]
cells_df['NFE2L2']=transpose_counts[adata.var.index.get_loc("NFE2L2")]
cells_df['UBE2D4']=transpose_counts[adata.var.index.get_loc("UBE2D4")]
cells_df['YPEL2']=transpose_counts[adata.var.index.get_loc("YPEL2")]

# genes dataframe
gene_table = adata.var
gene_table["gene_ids"]=gene_table.index

varm_np = np.array(adata.varm["PCs"])
gene_table["PCs_1"] = varm_np[:, 0]
gene_table["PCs_2"] = varm_np[:, 1]

# Cells charts

# creating a row chart based on the Row Chart implementation to show the leiden clustering
row_chart = RowChart(
    title="leiden",
    param="leiden",
    position=[380, 270],
    size=[260, 260]
)

# configuring the row chart
row_chart.set_axis_properties("x", {"textSize": 13, "label": "", "tickfont": 10})


# creating a dot plot based on the DotPlot implementation to show the gene expression of selected gene markers
dot_plot = DotPlot(
    title="leiden",
    params=["leiden", "ARVCF", "DOK3", "FAM210B", "GBGT1", "NFE2L2", "UBE2D4", "YPEL2"],
    size=[400, 250],
    position=[10, 10]
)

# configuring the dot plot
dot_plot.set_axis_properties("x", {"label": "", "textSize": 13, "tickfont": 10})
dot_plot.set_axis_properties("y", {"label": "", "textSize": 13, "tickfont": 10})
dot_plot.set_axis_properties("ry", {"label": "", "textSize": 13, "tickfont": 10})
dot_plot.set_color_scale(log_scale=False)
dot_plot.set_color_legend(True, [40, 10])
dot_plot.set_fraction_legend(True, [0, 0])


# creating a histogram plot based on the HistogramPlot implementation to show the distribution of the number of genes per counts
histogram_plot = HistogramPlot(
    title="n_genes_by_counts",
    param="n_genes_by_counts",
    bin_number=17,
    display_min=0,
    display_max=2500,
    size=[360, 250],
    position=[10, 280]
)

# configuring the histogram
histogram_plot.set_x_axis(size=30, label="n_genes_by_counts", textsize=13, tickfont=10)
histogram_plot.set_y_axis(size=45, label="frequency", textsize=13, tickfont=10, rotate_labels=False)

# creating a scatter plot based on the ScatterPlot3D implementation to show the 3 PCA clustering components
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
)

# configuring the scatter plot
scatter_plot.set_color_by("leiden")

# Genes charts

# creating a scatter plot using the ScatterPlot implementation to show the 2 PCA clustering components
scatter2D_plot = ScatterPlot(
    title="PCs_1 x PCs_2",
    params=["PCs_1", "PCs_2"],
    size=[250, 250],
    position=[270, 10]
)

# configuring the scatter plot
scatter2D_plot.set_default_color("#377eb8")
scatter2D_plot.set_brush("poly")
scatter2D_plot.set_opacity(0.8)
scatter2D_plot.set_radius(2)
scatter2D_plot.set_color_legend(display=True, position=[20, 1])
scatter2D_plot.set_axis_properties("x", {"label": "PCs_1", "textSize": 13, "tickfont": 10})
scatter2D_plot.set_axis_properties("y", {"label": "PCs_2", "textSize": 13, "tickfont": 10})
scatter2D_plot.set_color_by("dispersions")

# creating a table using the TablePlot implementation to show all the data associated with the genes dataframe
table_plot = TablePlot(
    title="Data table",
    params= gene_table.columns.values.tolist() ,
    size=[250, 500],
    position=[10, 10]
)

# creating a histogram using the HistogramPlot implementation to show a histogram of the number of genes per cell distribution
histogram_plot_2 = HistogramPlot(
    title="n_cells",
    param="n_cells",
    bin_number=17,
    display_min=0,
    display_max=2500,
    size=[240, 240],
    position=[270, 280]
)

# configuring the histogram
histogram_plot_2.set_x_axis(size=30, label="n_cells", textsize=13, tickfont=10)
histogram_plot_2.set_y_axis(size=45, label="frequency", textsize=13, tickfont=10, rotate_labels=False)

# setting up and serving the MDV project
base = os.path.expanduser('~/mdv')
project_path = os.path.join(base, 'pbmc3k') # defining the location where the project metadata will be stored
p = MDVProject(os.path.expanduser(project_path), delete_existing=True)

# adding the two data sources to the project
p.add_datasource("cells", cells_df)
p.add_datasource("genes", gene_table)

# converting the chart implementation outputs to JSON and setting up the project view
list_charts_cells = []
list_charts_genes = []

# cells panel
list_charts_cells.append(row_chart.plot_data)
list_charts_cells.append(dot_plot.plot_data)
list_charts_cells.append(histogram_plot.plot_data)
list_charts_cells.append(scatter_plot.plot_data)

# genes panel
list_charts_genes.append(table_plot.plot_data)
list_charts_genes.append(scatter2D_plot.plot_data)
list_charts_genes.append(histogram_plot_2.plot_data)

# setting the config combining the two panels
view_config = {'initialCharts': {"cells": list_charts_cells, "genes":list_charts_genes}}

# creating the link between the two datasets so that selecting a subset of genes to add the expression in cells is enabled
p.add_rows_as_columns_link("cells","genes","gene_ids","Gene Expression")
p.add_rows_as_columns_subgroup("cells","genes","Gene scores",adata.X) #add the gene expression

# creating the link between the two datasets so that selecting a subset of genes to add the expression in cells is enabled
p.add_rows_as_columns_link("genes","cells","cell_id","Gene Expression per cell")
p.add_rows_as_columns_subgroup("genes","cells","Gene scores",adata.X.T) #add the gene expression 

# adding the view to the project configuration
p.set_view("default", view_config)
p.set_editable(True)

# serving the project
p.serve()