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
    title = "Example title"
    param = "param1" #param1 should be a numerical variable. When the variable is named as "param", it can only take one data field. If it was "params" it would take more than one.
    bin_number = 40
    display_min = float(data_frame[param].min()) # Convert to Python float to ensure compatibility with JSON serialization
    display_max = float(data_frame[param].max()) # Convert to Python float to ensure compatibility with JSON serialization
    size = [500, 300]
    position = [100, 10]
    
    x_axis_settings = {
        'size': 40,
        'label': "label1",
        'textsize': 11,
        'tickfont': 11
    }
    
    y_axis_settings = {
        'size': 35,
        'label': "frequency",
        'textsize': 11,
        'tickfont': 14,
        'rotate_labels': False
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
