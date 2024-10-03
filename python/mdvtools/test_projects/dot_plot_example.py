import os
import json
import pandas as pd
from mdvtools.mdvproject import MDVProject
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
    data_path = '/Users/mariak/Documents/MDVmk/MDV/python/mdvtools/data/data_cells.csv'
    view_name = "default"
    
    # Create project
    project = MDVProject(project_path, delete_existing=True)
    
    # Load data
    data_frame = load_data(data_path)

    # Set the correct data type to the "leiden" data source (imports as integer but it should be str to appear as a category)
    data_frame['leiden'] = data_frame['leiden'].apply(str)
    
    # Add datasource
    project.add_datasource(data_path, data_frame)
    
    # DotPlot parameters
    title = "Dot Plot Example"
    params = ["leiden", "ARVCF", "DOK3", "FAM210B", "GBGT1", "NFE2L2", "UBE2D4", "YPEL2"]
    size = [792, 472]
    position = [10, 10]

    colorscale = {
        'log': False
    }
    
    # Create plot
    dot_plot = create_dot_plot(title, params, size, position, colorscale)
    
    # Convert plot to JSON and set view
    dot_plot_json = convert_plot_to_json(dot_plot)
    dotplot_view = {'initialCharts': {data_path: [dot_plot_json]}}
    
    project.set_view(view_name, dotplot_view)
    project.set_editable(True)
    project.serve()

if __name__ == "__main__":
    main()