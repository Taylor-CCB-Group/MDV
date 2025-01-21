import os
import pandas as pd
import scanpy as sc
import  numpy as np
from mdvtools.mdvproject import MDVProject
from mdvtools.charts.scatter_plot import ScatterPlot
import json 


def load_data(path):
    #Load data from the specified CSV file.
    return pd.read_csv(path, low_memory=False)

def convert_plot_to_json(plot):
    #Convert plot data to JSON format.
    return json.loads(json.dumps(plot.plot_data, indent=2).replace("\\", ""))
    

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

    # Add PCA data to the dataframe
    pca_np = np.array(adata.obsm["X_pca"])
    cells_df["PCA 1"] = pca_np[:, 0]
    cells_df["PCA 2"] = pca_np[:, 1]

    # Update datasource with the new columns provided through the scanpy object
    project.set_column(datasource_name, "PCA 1", cells_df["PCA 1"])
    project.set_column(datasource_name, "PCA 2", cells_df["PCA 2"])

    # ScatterPlot parameters
    scatter_title = "PCA Scatter Plot Colored by the leiden clusters"
    scatter_params = ["PCA 1", "PCA 2"]
    scatter_size = [792, 472]
    scatter_position = [820, 10]
    scatter_color = 'leiden'  # Color by the leiden clusters
    scatter_x_axis_settings = {'size': 30, 'label': "PCA 1", 'textsize': 13, 'tickfont': 10}
    scatter_y_axis_settings = {'size': 45, 'label': "PCA 2", 'textsize': 13, 'tickfont': 10, 'rotate_labels': False}
    
    # Create scatter plot
    scatter_plot = create_scatter_plot(
        scatter_title, scatter_params, scatter_size, scatter_position, scatter_color, scatter_x_axis_settings, scatter_y_axis_settings
    )
    
    # Convert plots to JSON and set view
    scatter_plot_json = convert_plot_to_json(scatter_plot)
    
    view_config = {'initialCharts': {datasource_name: [scatter_plot_json]}}
    
    project.set_view(view_name, view_config)
    project.set_editable(True)
    project.serve()

if __name__ == "__main__":
    main()