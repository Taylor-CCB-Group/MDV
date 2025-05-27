import os
import json
import pandas as pd
from mdvtools.mdvproject import MDVProject
from mdvtools.charts.histogram_plot import HistogramPlot

def create_histogram_plot(title, param, bin_number, display_min, display_max, size, position, x_axis_settings, y_axis_settings):
    """Create and configure a HistogramPlot instance with the given parameters."""
    plot = HistogramPlot(
        title=title,
        param=param, # When the variable is named as "param", it can only take one data field. If it was "params" it would take more than one.
        bin_number=bin_number,
        display_min=display_min,
        display_max=display_max,
        size=size,
        position=position
    )
    
    plot.set_x_axis(**x_axis_settings)  # x-axis settings
    plot.set_y_axis(**y_axis_settings)  # y-axis settings
    
    return plot

def load_data(path):
    """Load data from the specified CSV file."""
    return pd.read_csv(path, low_memory=False)

def convert_plot_to_json(plot):
    """Convert plot data to JSON format."""
    return json.loads(json.dumps(plot.plot_data, indent=2).replace("\\", ""))

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
    
    # Add datasource
    project.add_datasource(datasource_name, data_frame)
    
    # HistogramPlot parameters
    title = "Data Distribution"
    param = "param2"  # Changing the parameter name to "param2"
    bin_number = 75
    display_min = float(data_frame[param].min())  # Display min value for param2
    display_max = float(data_frame[param].max())  # Display max value for param2
    size = [1024, 512]
    position = [15, 20]

    x_axis_settings = {
        'size': 35,
        'label': "Value of param2",
        'textsize': 15,
        'tickfont': 12
    }

    y_axis_settings = {
        'size': 50,
        'label': "Frequency of param2",
        'textsize': 14,
        'tickfont': 12,
        'rotate_labels': True
    }

    
    # Create and configure plot
    histogram_plot = create_histogram_plot(
        title, param, bin_number, display_min, display_max, size, position, x_axis_settings, y_axis_settings
    )
    
    # Convert plot to JSON and set view
    histogram_chart_json = convert_plot_to_json(histogram_plot)
    histogram_view = {'initialCharts': {datasource_name: [histogram_chart_json]}}
    
    project.set_view(view_name, histogram_view)
    project.set_editable(True)
    project.serve()

if __name__ == "__main__":
    main()
