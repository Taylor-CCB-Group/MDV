import os
import pandas as pd
import scanpy as sc
import  numpy as np
from mdvtools.mdvproject import MDVProject
from mdvtools.charts.scatter_plot import ScatterPlot
import json 

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
    view_name = "Scatter Plot: Inflammation Score vs Age"
    
    # Load data
    data_path = "file_path"
    adata = sc.read_h5ad(data_path)
    cells_df = pd.DataFrame(adata.obs)
    cells_df.name = 'cells'
    
    # Add UMAP data to the dataframe
    umap_np = np.array(adata.obsm["X_umap"])
    cells_df["UMAP 1"] = umap_np[:, 0]
    cells_df["UMAP 2"] = umap_np[:, 1]
    
    # Create project
    project = MDVProject(project_path, delete_existing=True)
    
    # Add datasource
    project.add_datasource('cells', cells_df)

    # Update datasource with the new columns provided through the scanpy object
    project.set_column('cells', "UMAP 1", cells_df["UMAP 1"])
    project.set_column('cells', "UMAP 2", cells_df["UMAP 2"])
    
    # ScatterPlot parameters
    scatter_title = "Scatter Plot: Inflammation Score vs Age"
    scatter_params = ["Inflammation_score", "Age"]
    scatter_size = [792, 472]
    scatter_position = [820, 10]
    scatter_color = 'final_analysis'  # Assuming 'final_analysis' differentiates clusters
    scatter_x_axis_settings = {'size': 30, 'label': "Inflammation Score", 'textsize': 13, 'tickfont': 10}
    scatter_y_axis_settings = {'size': 45, 'label': "Age", 'textsize': 13, 'tickfont': 10, 'rotate_labels': False}
    
    # Create scatter plot
    scatter_plot = create_scatter_plot(
        scatter_title, scatter_params, scatter_size, scatter_position, scatter_color, scatter_x_axis_settings, scatter_y_axis_settings
    )
    
    # Convert plot to JSON and set view
    scatter_plot_json = convert_plot_to_json(scatter_plot)
    
    view_config = {'initialCharts': {'cells': [scatter_plot_json]}}
    
    project.set_view(view_name, view_config)
    project.set_editable(True)
    project.serve()

if __name__ == "__main__":
    main()