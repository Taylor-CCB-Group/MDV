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
    data_path = 'path_to_data'
    view_name = "default"
    datasource_name = "datasource_name"
    
    # Create project
    project = MDVProject(project_path, delete_existing=True)
    
    # Load data
    data_frame = load_data(data_path)

    # Add datasource
    project.add_datasource(datasource_name, data_frame)

    # TablePlot parameters
    title="Example title",
    params = ["param1", "param2", "param3", "param4", "param5", "param6"]  
    size = [850, 550]  
    position = [40, 40]  

    # Create plot
    table_plot = create_table_plot(title, params, size, position)
    
    # Convert plot to JSON and set view
    table_plot_json = convert_plot_to_json(table_plot)
    tableplot_view = {'initialCharts': {datasource_name: [table_plot_json]}}
    
    project.set_view(view_name, tableplot_view)
    project.set_editable(True)
    project.serve()

if __name__ == "__main__":
    main()
