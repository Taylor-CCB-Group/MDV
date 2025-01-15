import os
import json
import pandas as pd
from mdvtools.mdvproject import MDVProject
from mdvtools.charts.row_chart import RowChart

def create_row_chart(title, param, size, position, axis_settings):
    """Create and configure a RowChart instance with the given parameters."""
    plot = RowChart(
        title=title,
        param=param, #the param has to be just one categorical variable in the form of a string
        size=size,
        position=position
    )
    
    plot.set_axis_properties("x", axis_settings)  # axis settings
    
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
    data_frame = load_data(data_path)

    # Set the correct data type to the "leiden" data source (imports as integer but it should be str to appear as a category)
    data_frame['leiden'] = data_frame['leiden'].apply(str)
    
    # Add datasource
    project.add_datasource(datasource_name, data_frame)
    
    # RowChart parameters
    title = "Row Chart Example"
    param = "leiden"
    size = [792, 472]
    position = [10, 10]
    
    axis_settings = {
        'textSize': 13,
        'label': "Axis label",
        'tickfont': 10
    }
    
    # Create plot
    row_chart = create_row_chart(title, param, size, position, axis_settings)

    # Convert plot to JSON and set view
    row_chart_json = convert_plot_to_json(row_chart)
    rowchart_view = {'initialCharts': {datasource_name: [row_chart_json]}}
    
    project.set_view(view_name, rowchart_view)
    project.set_editable(True)
    project.serve()

if __name__ == "__main__":
    main()
