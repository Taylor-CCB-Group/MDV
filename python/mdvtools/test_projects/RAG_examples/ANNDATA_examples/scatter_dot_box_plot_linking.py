import os
import pandas as pd
import scanpy as sc
import  numpy as np
from mdvtools.mdvproject import MDVProject
from mdvtools.charts.dot_plot import DotPlot
from mdvtools.charts.box_plot import BoxPlot
from mdvtools.charts.scatter_plot import ScatterPlot
import json


def create_scatter_plot(title, params, size, position, color, brush, opacity, radius, legend_display, legend_position, xaxis_properties, yaxis_properties):
    """Create and configure a ScatterPlot instance with the given parameters."""
    plot = ScatterPlot(
        title=title,
        params=params,
        size=size,
        position=position
    )

    plot.set_default_color(color)
    plot.set_brush(brush)
    plot.set_opacity(opacity)
    plot.set_radius(radius)
    plot.set_color_legend(legend_display, legend_position)
    plot.set_axis_properties("x", xaxis_properties)
    plot.set_axis_properties("y", yaxis_properties)
    
    return plot

def create_dot_plot(title, params, size, position, colorscale):
    """Create and configure a DotPlot instance with the given parameters."""
    plot = DotPlot(
        title=title,
        params=params,
        size=size,
        position=position
    )
    plot.set_color_scale(colorscale)
    
    return plot

def create_box_plot(title, params, size, position):
    """Create and configure a BoxPlot instance with the given parameters."""
    plot = BoxPlot(
        title=title,
        params=params,
        size=size,
        position=position
    )
    return plot

def convert_plot_to_json(plot):
    return json.loads(json.dumps(plot.plot_data, indent=2).replace("\\", ""))

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
    scatter_title = "Scatter Plot"
    scatter_params = ["param3", "param4"]
    scatter_size = [792, 472]
    scatter_position = [820, 10]

    scatter_color = "param6"
    scatter_brush = "poly"
    scatter_opacity = 0.8
    scatter_radius = 0.2

    scatter_legend_display = True
    scatter_legend_position = [375,1]
              
    
    scatter_xaxis_properties = {"label": "label 1", 
             "textSize": 13, 
             "tickfont": 10
    }

    scatter_yaxis_properties = {"label": "label 2", 
             "textSize": 13, 
             "tickfont": 10
    }
    
    # Create scatter plot
    scatter_plot = create_scatter_plot(scatter_title, scatter_params, scatter_size, scatter_position, scatter_color, scatter_brush, scatter_opacity, scatter_radius, scatter_legend_display, scatter_legend_position, scatter_xaxis_properties, scatter_yaxis_properties)

    
    # DotPlot parameters
    param2 = "param2"
    param2_index = data_frame_var['name'].tolist().index(param2)

    # DotPlot parameters
    dot_title = "Dot plot example title"
    dot_params = ["param1", f"gs|{param2}(gs)|{param2_index}"]
    dot_size = [1000, 250]
    dot_position = [10, 500]

    dot_colorscale = {'log': False}
    
    # Create dot plot
    dot_plot = create_dot_plot(dot_title, dot_params, dot_size, dot_position, dot_colorscale)
    
    # BoxPlot parameters
    param2 = "param2"
    param2_index = data_frame_var['name'].tolist().index(param2)
    box_title = "Example title"
    box_params = ["param4", f"gs|{param2}(gs)|{param2_index}"]
    box_size = [1015, 557]
    box_position = [500, 500]
    
    # Create box plot
    box_plot = create_box_plot(
        box_title, box_params, box_size, box_position)
    
    # Convert plots to JSON and set view
    scatter_plot_json = convert_plot_to_json(scatter_plot)
    dot_plot_json = convert_plot_to_json(dot_plot)
    box_plot_json = convert_plot_to_json(box_plot)
    
    view_config = {'initialCharts': {datasource_name: [scatter_plot_json, dot_plot_json, box_plot_json]}}

    project.set_view(view_name, view_config)
    project.set_editable(True)
    project.serve()

if __name__ == "__main__":
    main()