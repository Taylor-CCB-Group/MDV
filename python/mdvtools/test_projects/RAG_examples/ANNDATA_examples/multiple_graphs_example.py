import os
import json
import pandas as pd
import scanpy as sc
from mdvtools.mdvproject import MDVProject
from mdvtools.charts.box_plot import BoxPlot

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
    
    # BoxPlot parameters
    title = "Box Plot Example"
    params_pairs = [["param1", "param2"], ["param1", "param3"], ["param1", "param3"]] 
    size = [350, 350]
    initial_position = [10, 10]

    list_boxplot_charts_json = []
    
    # Create the multiple plots, convert them all to the json format and add them all to a list for integration to MDV
    for i in range(0, len(params_pairs)):
        plot = create_box_plot(title, params_pairs[i], size, [x + y for x, y in zip(initial_position, [k * i for k in [0, 360]])])
    
        # Convert plot to JSON and set view
        boxplot_charts_json = convert_plot_to_json(plot)
        list_boxplot_charts_json.append(boxplot_charts_json)
    
    boxplot_view = {'initialCharts': {datasource_name: list_boxplot_charts_json}}
    
    project.set_view(view_name, boxplot_view)
    project.set_editable(True)
    project.serve()

if __name__ == "__main__":
    main()