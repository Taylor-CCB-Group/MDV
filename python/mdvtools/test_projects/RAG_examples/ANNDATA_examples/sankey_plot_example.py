import os
import json
import pandas as pd
import scanpy as sc
from mdvtools.mdvproject import MDVProject
from mdvtools.charts.sankey_plot import SankeyPlot

def create_sankey_plot(title, params, size, position):
    """Create and configure a SankeyPlot instance with the given parameters."""
    plot = SankeyPlot(
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
    adata = sc.read_h5ad(data_path, backed='r')
    data_frame_obs = pd.DataFrame(adata.obs)

    data_frame_var = pd.DataFrame(adata.var)
    data_frame_var['name'] = adata.var_names.to_list()
    
    # Add datasource
    project.add_datasource(datasource_name, data_frame_obs)
    project.add_datasource(datasource_name_2, data_frame_var)

    # SankeyPlot parameters
    title="Example title",
    params = ["param1", "param2"] # The params should be two categorical variables.
    size = [792, 472]
    position = [10, 10]

    # Create plot
    sankey_plot = create_sankey_plot(title, params, size, position)
    
    # Convert plot to JSON and set view
    sankey_plot_json = convert_plot_to_json(sankey_plot)
    sankeyplot_view = {'initialCharts': {datasource_name: [sankey_plot_json]}}
    
    project.set_view(view_name, sankeyplot_view)
    project.set_editable(True)
    project.serve()

if __name__ == "__main__":
    main()
