import os
import json
import pandas as pd
from mdvtools.mdvproject import MDVProject
from mdvtools.charts.scatter_plot_3D import ScatterPlot3D

def create_scatter_plot(title, params, size, position, color, brush, opacity, radius, camera, center, category_color):
    """Create and configure a ScatterPlot3D instance with the given parameters."""
    plot = ScatterPlot3D(
        title=title,
        params=params,
        size=size,
        position=position
    )

    plot.set_default_color(color)
    plot.set_brush(brush)
    plot.set_opacity(opacity)
    plot.set_radius(radius)
    plot.set_camera(camera)
    plot.set_center(center)
    plot.set_color_by(category_color)
    
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
    
    # ScatterPlot3D parameters
    title = "Scatter Plot Example"
    params = ["param1", "param2", "param3"] #param1, param2 and param3 should all be a numerical variables
    size = [700, 500]
    position = [50, 50]
    center = [0, 1950, 2050]

    color = "#9467bd"
    brush = "default"
    radius = 8.0
    opacity = 0.9
    camera = {"distance": 0.05, "theta": 0.08, "phi": 0.28}

    category_color = "categorial variable"

    # Create plot
    scatter_plot = create_scatter_plot(title, params, size, position, color, brush, opacity, radius, camera, center, category_color)
    
    # Convert plot to JSON and set view
    scatter_plot_json = convert_plot_to_json(scatter_plot)
    scatterplot3D_view = {'initialCharts': {datasource_name: [scatter_plot_json]}}
    
    project.set_view(view_name, scatterplot3D_view)
    project.set_editable(True)
    project.serve()

if __name__ == "__main__":
    main()
