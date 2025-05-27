import os
import json
import pandas as pd
import scanpy as sc
from mdvtools.mdvproject import MDVProject
from mdvtools.conversions import convert_scanpy_to_mdv
from mdvtools.charts.dot_plot import DotPlot

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
    datasource_name_2 = "datasource_name_2"
    
    # Load data
    adata = sc.read_h5ad(data_path)
    data_frame_obs = pd.DataFrame(adata.obs)

    data_frame_var = pd.DataFrame(adata.var)
    data_frame_var['name'] = adata.var_names.to_list()

    # Create project
    project = convert_scanpy_to_mdv(project_path, adata)
        
    # Add datasource
    project.add_datasource(datasource_name, data_frame_obs)
    project.add_datasource(datasource_name_2, data_frame_var)
    
    # DotPlot parameters
    title = "Example Title"
    params = ["param1", "param2", "param3", "param4", "param5", "param6", "param7", "param8"] # The 'params' list can accept any number of arguments, param1 should be a categorical variable and the rest of params should be numerical variables.
    size = [792, 472]
    position = [10, 10]

    colorscale = {
        'log': False
    }
    
    # Create plot
    dot_plot = create_dot_plot(title, params, size, position, colorscale)
    
    # Convert plot to JSON and set view
    dot_plot_json = convert_plot_to_json(dot_plot)
    dotplot_view = {'initialCharts': {datasource_name: [dot_plot_json]}}
    
    project.set_view(view_name, dotplot_view)
    project.set_editable(True)
    project.serve()

if __name__ == "__main__":
    main()