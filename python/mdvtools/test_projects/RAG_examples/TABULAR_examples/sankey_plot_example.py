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
    data_frame = pd.DataFrame(data_path)

    # Add datasource
    project.add_datasource(datasource_name, data_frame)

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
