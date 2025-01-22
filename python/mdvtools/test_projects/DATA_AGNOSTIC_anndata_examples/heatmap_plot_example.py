import os
import json
import pandas as pd
import scanpy as sc
from mdvtools.mdvproject import MDVProject
from mdvtools.charts.heatmap_plot import HeatmapPlot

def create_heatmap_plot(title, params, size, position, colorscale, x_axis_settings, y_axis_settings):
    """Create and configure a HeatmapPlot instance with the given parameters."""
    plot = HeatmapPlot(
        title=title,
        params=params,
        size=size,
        position=position
    )
    
    plot.set_color_scale(colorscale)  # colorscale setting
    plot.set_x_axis(**x_axis_settings)  # x-axis settings
    plot.set_y_axis(**y_axis_settings)  # y-axis settings
    
    return plot

def convert_plot_to_json(plot):
    """Convert plot data to JSON format."""
    return json.loads(json.dumps(plot.plot_data, indent=2).replace("\\\\", ""))

def main():
    """Main function to create the project and serve it."""
    # Constants
    project_path = os.path.expanduser('~/mdv/project')
    data_path = "path_to_data"
    view_name = "default"
    datasource_name = "datasource_name"
    
    # Create project
    project = MDVProject(project_path, delete_existing=True)
    
    # Load data
    adata = sc.read_h5ad(data_path)
    data_frame = pd.DataFrame(adata.obs)
    
    # Add datasource
    project.add_datasource(datasource_name, data_frame)
    
    # HeatmapPlot parameters
    title = "Example title"
    params = ["param1", "param2", "param3", "param4", "param5", "param6", "param7", "param8"] # The 'params' list can accept any number of arguments, param1 should be a categorical variable and the rest of params should be numerical variables.
    size = [792, 472]
    position = [10, 10]
    
    colorscale = {
        'log': False
    }
    
    x_axis_settings = {
        'axis_labels': "X-Axis Label",
        'axis_title': "X-Axis Title",
    }
    
    y_axis_settings = {
        'axis_labels': "Y-Axis Label",
        'axis_title': "Y-Axis Title",
        #'rotate_labels': 45
    }
    
    # Create plot
    heatmap_plot = create_heatmap_plot(title, params, size, position, colorscale, x_axis_settings, y_axis_settings)

    # Convert plot to JSON and set view
    heatmap_plot_json = convert_plot_to_json(heatmap_plot)
    heatmap_view = {'initialCharts': {datasource_name: [heatmap_plot_json]}}
    
    project.set_view(view_name, heatmap_view)
    project.set_editable(True)
    project.serve()

if __name__ == "__main__":
    main()
