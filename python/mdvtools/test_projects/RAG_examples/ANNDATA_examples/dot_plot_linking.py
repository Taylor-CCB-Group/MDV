import os
import pandas as pd
import scanpy as sc
from mdvtools.mdvproject import MDVProject
from mdvtools.charts.dot_plot import DotPlot
import json 


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

    # DotPlot parameters
    param2 = "param2"
    param2_index = data_frame_var['name'].tolist().index(param2)

    # DotPlot parameters
    dot_title = "Dot plot example title"
    dot_params = ["param1", f"gs|{param2}(gs)|{param2_index}"]
    dot_size = [800, 250]
    dot_position = [10, 500]

    dot_colorscale = {'log': False}
    
    # Create dot plot
    dot_plot = create_dot_plot(dot_title, dot_params, dot_size, dot_position, dot_colorscale)

    dot_plot_json = convert_plot_to_json(dot_plot)
    
    view_config = {'initialCharts': {datasource_name: [dot_plot_json]}}

    project.set_view(view_name, view_config)
    project.set_editable(True)
    project.serve()

if __name__ == "__main__":
    main()