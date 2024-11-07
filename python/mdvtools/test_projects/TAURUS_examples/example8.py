## This file creates a new MDV project named 'project' containing a view named 'default'
## The datasource used is an h5ad file that was provided at run time
## The view 'default' shows one stacked bar plot, one dot plot, one scatter plot and one box plot.
## The h5ad file is an AnnData object and here the obs attribute was used.
## The stacked bar plot shows the distribution of various cell states across different patients.
## By examining the cell state proportions for each patient, one can observe if there are unique patterns or 
## distributions specific to certain patients. For example, patients with particular conditions may show a higher abundance of certain cell states.
## UMAP Scatter Plot uses UMAP 1 and UMAP 2 as coordinates. Points are colored by cell state to show different cell types or states.
## This UMAP plot reveals clusters of cells based on similarity, with each cluster representing specific cell types or states. 
## It helps visualize the cellular landscape and detect distinct cell populations.
## The dot plot uses cell state on the x-axis, with n_genes_by_counts and total_counts as variables represented by dot size and color intensity.
## The dot plot gives an overview of gene expression levels per cell state. It helps identify cell types with high or low expression levels, indicating variations in cellular activity across states.
## The box plot shows Disease on the x-axis and n_genes_by_counts on the y-axis.
## The box plot compares gene expression levels across different disease states, highlighting differences in gene expression that may relate to disease pathology or severity.

import os
import pandas as pd
import scanpy as sc
import sys
import  numpy as np
from mdvtools.mdvproject import MDVProject
from mdvtools.charts.dot_plot import DotPlot
from mdvtools.charts.box_plot import BoxPlot
from mdvtools.charts.scatter_plot import ScatterPlot
from mdvtools.charts.stacked_row_plot import StackedRowChart
import json 

def create_stacked_row_plot(title, params, size, position, legend_display, xaxis_properties, yaxis_properties):
    plot = StackedRowChart(
        title=title,
        params=params,
        size=size,
        position=position
    )
    plot.set_color_legend(legend_display)
    plot.set_axis_properties("x", xaxis_properties)
    plot.set_axis_properties("y", yaxis_properties)
    return plot

def create_scatter_plot(title, params, size, position, color, x_axis_settings, y_axis_settings):
    plot = ScatterPlot(
        title=title,
        params=params,
        size=size,
        position=position
    )
    plot.set_color_by(color)
    plot.set_axis_properties("x", x_axis_settings)
    plot.set_axis_properties("y", y_axis_settings)
    return plot

def create_dot_plot(title, params, size, position):
    plot = DotPlot(
        title=title,
        params=params,
        size=size,
        position=position
    )
    plot.set_axis_properties("x", {"label": "", "textSize": 13, "tickfont": 10})
    plot.set_axis_properties("y", {"label": "", "textSize": 13, "tickfont": 10})
    plot.set_color_scale(log_scale=False)
    plot.set_color_legend(True, [40, 10])
    plot.set_fraction_legend(True, [0, 0])
    return plot

def create_box_plot(title, params, size, position, plot_id):
    plot = BoxPlot(
        title=title,
        params=params,
        size=size,
        position=position,
        id=plot_id
    )

    return plot

def convert_plot_to_json(plot):
    return json.loads(json.dumps(plot.plot_data, indent=2).replace("\\", ""))

def main():
    project_path = os.path.expanduser('~/mdv/project')
    view_name = "default"
    
    # Load data
    adata = sc.read_h5ad(sys.argv[1])
    cells_df = pd.DataFrame(adata.obs)
    cells_df.name = 'cells'
    
    # Rename 'final_analysis' to 'cell state'
    cells_df.rename(columns={"final_analysis": "cell state"}, inplace=True)
    
    # Add UMAP data to the dataframe
    umap_np = np.array(adata.obsm["X_umap"])
    cells_df["UMAP 1"] = umap_np[:, 0]
    cells_df["UMAP 2"] = umap_np[:, 1]
    
    # Create project
    project = MDVProject(project_path, delete_existing=True)
    
    # Add datasource
    project.add_datasource('cells', cells_df)
    
    # StackedRowChart parameters
    stacked_title = "Abundance of Cell States per Patient"
    stacked_params = ["Patient", "cell state"]
    stacked_size = [792, 472]
    stacked_position = [10, 10]
    stacked_legend_display = True
    stacked_xaxis_properties = {"label": "Patient", "textSize": 13, "tickfont": 10}
    stacked_yaxis_properties = {"label": "Cell State", "textSize": 13, "tickfont": 10}
    
    # Create stacked row plot
    stacked_row_plot = create_stacked_row_plot(
        stacked_title, stacked_params, stacked_size, stacked_position, stacked_legend_display, stacked_xaxis_properties, stacked_yaxis_properties
    )
    
    # ScatterPlot parameters
    scatter_title = "UMAP Scatter Plot"
    scatter_params = ["UMAP 1", "UMAP 2"]
    scatter_size = [792, 472]
    scatter_position = [820, 10]
    scatter_color = 'cell state'
    scatter_x_axis_settings = {'size': 30, 'label': "UMAP 1", 'textsize': 13, 'tickfont': 10}
    scatter_y_axis_settings = {'size': 45, 'label': "UMAP 2", 'textsize': 13, 'tickfont': 10, 'rotate_labels': False}
    
    # Create scatter plot
    scatter_plot = create_scatter_plot(
        scatter_title, scatter_params, scatter_size, scatter_position, scatter_color, scatter_x_axis_settings, scatter_y_axis_settings
    )
    
    # DotPlot parameters
    dot_title = "Gene expression counts per cell state dot plot"
    dot_params = ["cell state", "n_genes_by_counts", "total_counts"]
    dot_size = [400, 250]
    dot_position = [10, 500]
    
    # Create dot plot
    dot_plot = create_dot_plot(dot_title, dot_params, dot_size, dot_position)
    
    # BoxPlot parameters
    box_title = "Gene expression counts per disease box plot"
    box_params = ["Disease", "n_genes_by_counts"]
    box_size = [615, 557]
    box_position = [500, 500]
    box_plot_id = "boxPlot1"
    
    # Create box plot
    box_plot = create_box_plot(
        box_title, box_params, box_size, box_position, box_plot_id)
    
    # Convert plots to JSON and set view
    stacked_row_plot_json = convert_plot_to_json(stacked_row_plot)
    scatter_plot_json = convert_plot_to_json(scatter_plot)
    dot_plot_json = convert_plot_to_json(dot_plot)
    box_plot_json = convert_plot_to_json(box_plot)
    
    view_config = {'initialCharts': {'cells': [stacked_row_plot_json, scatter_plot_json, dot_plot_json, box_plot_json]}}
    
    project.set_view(view_name, view_config)
    project.set_editable(True)
    project.serve()

if __name__ == "__main__":
    main()