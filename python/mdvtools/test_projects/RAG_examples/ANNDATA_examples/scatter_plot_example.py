import os
import json
import pandas as pd
import scanpy as sc
from mdvtools.mdvproject import MDVProject
from mdvtools.charts.scatter_plot import ScatterPlot

def create_scatter_plot(title, params, size, position, color, brush, opacity, radius, legend_display, legend_position, xaxis_properties, yaxis_properties):
    """Create and configure a ScatterPlot instance with the given parameters."""
    plot = ScatterPlot(
        title=title,
        params=params,
        size=size,
        position=position
    )

    plot.set_default_color(color)
    plot.set_color_by(color)
    plot.set_brush(brush)
    plot.set_opacity(opacity)
    plot.set_radius(radius)
    plot.set_color_legend(legend_display, legend_position)
    plot.set_axis_properties("x", xaxis_properties)
    plot.set_axis_properties("y", yaxis_properties)
    
    return plot

def convert_plot_to_json(plot):
    """Convert plot data to JSON format."""
    return json.loads(json.dumps(plot.plot_data, indent=2).replace("\\\\", ""))

def main():
    """Main function to create the project and serve it."""
    # Constants
    project_path = os.path.expanduser('~/mdv/project')
    data_path = "path_to_data"
    view_name = "Example view name"
    datasource_name = "datasource_name"
    datasource_name_2 = "datasource_name_2"
    
    # Create project
    project = MDVProject(project_path, delete_existing=True)
    
    # Load data
    adata = sc.read_h5ad(data_path, backed='r')
    data_frame_obs = pd.DataFrame(adata.obs)

    data_frame_var = pd.DataFrame(adata.var)
    data_frame_var['name'] = adata.var_names.to_list()
    
    # Add datasource
    project.add_datasource(datasource_name, data_frame_obs)
    project.add_datasource(datasource_name_2, data_frame_var)
    
    # ScatterPlot parameters
    title = "Scatter Plot Example"
    params = ["param1", "param2"] #param1 and param2 should all be a numerical variables
    size = [792, 472]
    position = [10, 10]

    color = "param3"
    brush = "poly"
    opacity = 0.8
    radius = 4

    legend_display = True
    legend_position = [375,1]
              
    
    xaxis_properties = {"label": "label 1", 
             "textSize": 13, 
             "tickfont": 10
    }

    yaxis_properties = {"label": "label 2", 
             "textSize": 13, 
             "tickfont": 10
    }

    # Create plot
    scatter_plot = create_scatter_plot(title, params, size, position, color, brush, opacity, radius, legend_display, legend_position, xaxis_properties, yaxis_properties)
    
    # Convert plot to JSON and set view
    scatter_plot_json = convert_plot_to_json(scatter_plot)
    scatterplot_view = {'initialCharts': {datasource_name: [scatter_plot_json]}}
    
    project.set_view(view_name, scatterplot_view)
    project.set_editable(True)
    project.serve()

if __name__ == "__main__":
    main()
