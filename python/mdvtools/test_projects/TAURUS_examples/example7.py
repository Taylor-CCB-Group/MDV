## This file creates a new MDV project named 'project' containing a view named 'default'
## The datasource used is an h5ad file that was provided at run time
## The view 'default' shows one stacked bar plot and one scatter plot
## The h5ad file is an AnnData object and here the obs and obsm attributes were used
## The stacked bar plot shows the distribution of various cell states across different patients.
## By examining the cell state proportions for each patient, one can observe if there are unique patterns or 
## distributions specific to certain patients. For example, patients with particular conditions may show a higher abundance of certain cell states.
## The UMAP scatterplot shows clusters of cells grouped by similarity, which can reveal distinct cell populations or subtypes.
## UMAP is a dimensionality reduction technique commonly used in single-cell analysis to visualize high-dimensional data.


import os
import pandas as pd
import scanpy as sc
import numpy as np
import sys
from mdvtools.mdvproject import MDVProject
from mdvtools.charts.scatter_plot import ScatterPlot
from mdvtools.charts.stacked_row_plot import StackedRowChart
import json 

    
def create_stacked_row_plot(title, params, size, position, legend_display, xaxis_properties, yaxis_properties):
    """Create and configure a StackedRowChart instance with the given parameters."""
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

def create_scatter_plot(title, params, size, position, x_axis_settings, y_axis_settings):
    """Create and configure a ScatterPlot instance with the given parameters."""
    plot = ScatterPlot(
        title=title,
        params=params,
        size=size,
        position=position
    )

    plot.set_axis_properties("x", x_axis_settings)
    plot.set_axis_properties("y", y_axis_settings)

    return plot

def convert_plot_to_json(plot):
    """Convert plot data to JSON format."""
    return json.loads(json.dumps(plot.plot_data, indent=2).replace("\\", ""))

def main():
    """Main function to create the project and serve it."""
    # Constants
    project_path = os.path.expanduser('~/mdv/project')
    view_name = "default"
    
    # Load data
    data_path = "file_path"
    adata = sc.read_h5ad(data_path)
    cells_df = pd.DataFrame(adata.obs)
    datasource_name = "datasource_name"
    
    # Add UMAP data to the dataframe
    umap_np = np.array(adata.obsm["X_umap"])
    cells_df["UMAP 1"] = umap_np[:, 0]
    cells_df["UMAP 2"] = umap_np[:, 1]
    
    # Create project
    project = MDVProject(project_path, delete_existing=True)
    
    # Add datasource
    project.add_datasource(datasource_name, cells_df)

    # Update datasource with the new columns provided through the scanpy object
    project.set_column(datasource_name, "UMAP 1", cells_df["UMAP 1"])
    project.set_column(datasource_name, "UMAP 2", cells_df["UMAP 2"])
    
    # StackedRowChart parameters
    stacked_title = "Abundance of Cell States per Patient"
    stacked_params = ["final_analysis", "Patient"]
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
    
    scatter_x_axis_settings = {'size': 30, 'label': "UMAP 1", 'textsize': 13, 'tickfont': 10}
    scatter_y_axis_settings = {'size': 45, 'label': "UMAP 2", 'textsize': 13, 'tickfont': 10, 'rotate_labels': False}
    
    # Create scatter plot
    scatter_plot = create_scatter_plot(
        scatter_title, scatter_params, scatter_size, scatter_position, scatter_x_axis_settings, scatter_y_axis_settings
    )
    
    # Convert plots to JSON and set view
    stacked_row_plot_json = convert_plot_to_json(stacked_row_plot)
    scatter_plot_json = convert_plot_to_json(scatter_plot)
    
    view_config = {'initialCharts': {datasource_name: [stacked_row_plot_json, scatter_plot_json]}}
    
    project.set_view(view_name, view_config)
    project.set_editable(True)
    project.serve()

if __name__ == "__main__":
    main()