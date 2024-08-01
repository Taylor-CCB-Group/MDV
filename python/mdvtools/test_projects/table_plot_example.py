import os
import json
import pandas as pd
from mdvtools.mdvproject import MDVProject
from mdvtools.charts.table_plot import TablePlot

def create_table_plot(title, params, size, position):
    """Create and configure a TablePlot instance with the given parameters."""
    plot = TablePlot(
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
    data_path = '/Users/mariak/Documents/MDVmk/MDV/python/mdvtools/llm/sample_data/data_cells.csv'
    view_name = "default"
    
    # Create project
    project = MDVProject(project_path, delete_existing=True)
    
    # Load data
    data_frame = load_data(data_path)

    # Add datasource
    project.add_datasource(data_path, data_frame)

    # TablePlot parameters
    title="sample_id",
    params = ["leiden", "ARVCF", "DOK3", "FAM210B", "GBGT1", "NFE2L2", "UBE2D4", "YPEL2"]
    size = [792, 472]
    position = [10, 10]

    # Create plot
    table_plot = create_table_plot(title, params, size, position)
    
    # Convert plot to JSON and set view
    table_plot_json = convert_plot_to_json(table_plot)
    tableplot_view = {'initialCharts': {data_path: [table_plot_json]}}
    
    project.set_view(view_name, tableplot_view)
    project.set_editable(True)
    project.serve()

if __name__ == "__main__":
    main()
