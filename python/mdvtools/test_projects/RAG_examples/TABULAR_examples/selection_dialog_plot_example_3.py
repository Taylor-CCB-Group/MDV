import os
import json
import pandas as pd
from mdvtools.mdvproject import MDVProject
from mdvtools.charts.selection_dialog_plot import SelectionDialogPlot

def create_selection_dialog_plot(title, params, size, position):
    """Create and configure a SelectionDialogPlot instance with the given parameters."""
    plot = SelectionDialogPlot(
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

    # SelectionDialogPlot parameters
    title="Example title",
    params = ["param1", "param2", "param3"] # The 'param' list can accept any number of arguments, where each argument can be of any type, including categorical or variable types.
    size = [896, 340]
    position = [100, 10]

    # Create plot
    selection_dialog_plot = create_selection_dialog_plot(title, params, size, position)
    
    # Convert plot to JSON and set view
    selection_dialog_plot_json = convert_plot_to_json(selection_dialog_plot)
    selectiondialogplot_view = {'initialCharts': {datasource_name: [selection_dialog_plot_json]}}
    
    project.set_view(view_name, selectiondialogplot_view)
    project.set_editable(True)
    project.serve()

if __name__ == "__main__":
    main()
