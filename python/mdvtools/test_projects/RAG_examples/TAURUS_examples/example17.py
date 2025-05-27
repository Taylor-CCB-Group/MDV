import os
import pandas as pd
import scanpy as sc
import  numpy as np
from mdvtools.mdvproject import MDVProject
from mdvtools.charts.scatter_plot import ScatterPlot
from mdvtools.charts.stacked_row_plot import StackedRowChart
import json 

def create_stacked_row_plot(title, params, size, position, color_legend, xaxis_properties, yaxis_properties):
    plot = StackedRowChart(
        title=title,
        params=params,
        size=size,
        position=position
    )
    plot.set_color_legend(color_legend["display"], color_legend["pos"])
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

def convert_plot_to_json(plot):
    return json.loads(json.dumps(plot.plot_data, indent=2).replace("\\", ""))

def main():
    project_path = os.path.expanduser('~/mdv/project')
    view_name = "default"
    
    # Load data
    data_path = "file_path"
    adata = sc.read_h5ad(data_path)
    cells_df = pd.DataFrame(adata.obs)
    datasource_name = "datasource_name"
    
    # Create project
    project = MDVProject(project_path, delete_existing=True)
    
    # Add datasource
    project.add_datasource(datasource_name, cells_df)

    # Add UMAP data to the dataframe
    umap_np = np.array(adata.obsm["X_umap"])
    cells_df["UMAP 1"] = umap_np[:, 0]
    cells_df["UMAP 2"] = umap_np[:, 1]

    # Update datasource with the new columns provided through the scanpy object
    project.set_column(datasource_name, "UMAP 1", cells_df["UMAP 1"])
    project.set_column(datasource_name, "UMAP 2", cells_df["UMAP 2"])

    # StackedRowChart parameters
    stacked_title = "Distribution of Cell States Across Conditions"
    stacked_params = ["Remission_status", "bucket"]
    stacked_size = [792, 472]
    stacked_position = [10, 10]
    stacked_color_legend = {"display" : True,
                    "pos" : [375,1]}
    stacked_xaxis_properties = {"label": "Remission Status", "textSize": 13, "tickfont": 10}
    stacked_yaxis_properties = {"label": "Cell State", "textSize": 13, "tickfont": 10}
    
    # Create stacked row plot
    stacked_row_plot = create_stacked_row_plot(
        stacked_title, stacked_params, stacked_size, stacked_position, stacked_color_legend, stacked_xaxis_properties, stacked_yaxis_properties
    )
    
    # ScatterPlot parameters
    scatter_title = "Cell-Cell Interaction Network"
    scatter_params = ["UMAP 1", "UMAP 2"]
    scatter_size = [792, 472]
    scatter_position = [820, 10]
    scatter_color = 'major'
    scatter_x_axis_settings = {'size': 30, 'label': "UMAP 1", 'textsize': 13, 'tickfont': 10}
    scatter_y_axis_settings = {'size': 45, 'label': "UMAP 2", 'textsize': 13, 'tickfont': 10, 'rotate_labels': False}
    
    # Create scatter plot
    scatter_plot = create_scatter_plot(
        scatter_title, scatter_params, scatter_size, scatter_position, scatter_color, scatter_x_axis_settings, scatter_y_axis_settings
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