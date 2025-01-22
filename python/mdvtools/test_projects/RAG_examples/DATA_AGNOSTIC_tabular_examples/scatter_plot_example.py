import os
import json
import pandas as pd
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
    plot.set_brush(brush)
    plot.set_opacity(opacity)
    plot.set_radius(radius)
    plot.set_color_legend(legend_display, legend_position)
    plot.set_axis_properties("x", xaxis_properties)
    plot.set_axis_properties("y", yaxis_properties)
    
    return plot

def load_data(path):
    """Load data from the specified CSV file."""
    return pd.read_csv(path, low_memory=False)

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
    data_frame = load_data(data_path)
 
    # Add datasource
    project.add_datasource(datasource_name, data_frame)
    
    # ScatterPlot parameters
    title = "Scatter Plot Example"
    params = ["param1", "param2"] #param1 and param2 should all be a numerical variables
    size = [792, 472]
    position = [10, 10]

    color = "#377eb8"
    brush = "poly"
    opacity = 0.8
    radius = 0.2

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
