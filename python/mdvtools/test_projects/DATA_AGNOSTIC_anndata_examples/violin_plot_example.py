import os
import json
import pandas as pd
import scanpy as sc
from mdvtools.mdvproject import MDVProject
from mdvtools.charts.violin_plot import ViolinPlot

def create_violin_plot(title, params, size, position):
    """Create and configure a ViolinPlot instance with the given parameters."""
    plot = ViolinPlot(
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
    view_name = "default"
    datasource_name = "datasource_name"
    
    # Create project
    project = MDVProject(project_path, delete_existing=True)
    
    # Load data
    adata = sc.read_h5ad(data_path)
    data_frame = pd.DataFrame(adata.obs)
    
    # Add datasource
    project.add_datasource(datasource_name, data_frame)
    
    # ViolinPlot parameters
    title = "Violin Plot Example"
    params = ["param1", "param2"] #param1 should be a categorical variable and param2 should be a numerical variable.
    size = [792, 472]
    position = [10, 10]
    
    # Create plot
    plot = create_violin_plot(title, params, size, position)
    
    # Convert plot to JSON and set view
    ViolinPlot_chart_json = convert_plot_to_json(plot)
    ViolinPlot_view = {'initialCharts': {datasource_name: [ViolinPlot_chart_json]}}
    
    project.set_view(view_name, ViolinPlot_view)
    project.set_editable(True)
    project.serve()

if __name__ == "__main__":
    main()
