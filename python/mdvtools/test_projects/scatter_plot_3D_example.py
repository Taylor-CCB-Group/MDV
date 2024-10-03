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
    
    # Create project
    project = MDVProject(project_path, delete_existing=True)
    
    # Load data
    data_frame = load_data(data_path)

    # Set the correct data type to the "leiden" data source (imports as integer but it should be str to appear as a category)
    data_frame['leiden'] = data_frame['leiden'].apply(str)
    
    # Add datasource
    project.add_datasource(data_path, data_frame)
    
    # ScatterPlot3D parameters
    title = "Scatter Plot Example"
    params = ["ARVCF", "DOK3", "YPEL2"]
    size = [300, 300]
    position = [10, 10]
    center = [0,0,0]

    color = "#377eb8"
    brush = "default"
    radius = 0.9
    opacity = 0.8
    category_color = "leiden"
    
    camera={"distance": 37, "theta": -1.038, "phi": 0.261}


    # Create plot
    scatter_plot = create_scatter_plot(title, params, size, position, color, brush, opacity, radius, camera, center, category_color)
    
    # Convert plot to JSON and set view
    scatter_plot_json = convert_plot_to_json(scatter_plot)
    scatterplot3D_view = {'initialCharts': {data_path: [scatter_plot_json]}}
    
    project.set_view(view_name, scatterplot3D_view)
    project.set_editable(True)
    project.serve()

if __name__ == "__main__":
    main()
