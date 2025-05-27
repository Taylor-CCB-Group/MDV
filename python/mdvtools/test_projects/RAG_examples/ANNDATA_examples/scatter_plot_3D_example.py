import os
import json
import pandas as pd
import scanpy as sc
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
    adata = sc.read_h5ad(data_path)
    data_frame_obs = pd.DataFrame(adata.obs)

    data_frame_var = pd.DataFrame(adata.var)
    data_frame_var['name'] = adata.var_names.to_list()
    
    # Add datasource
    project.add_datasource(datasource_name, data_frame_obs)
    project.add_datasource(datasource_name_2, data_frame_var)
    
    # ScatterPlot3D parameters
    title = "Scatter Plot Example"
    params = ["param1", "param2", "param3"] #param1, param2 and param3 should all be a numerical variables
    size = [300, 300]
    position = [10, 10]
    center = [44.495425,-7.910625500000002,10.481393999999998]

    color = "param4"
    brush = "default"
    radius = 9.0
    opacity = 0.8
    category_color = "categorial variable"
    
    camera={"distance":660.6206799999997,"theta":-1.038,"phi":0.261}

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
