import os
import pandas as pd
from mdvtools.mdvproject import MDVProject
from mdvtools.charts.heatmap_plot import HeatmapPlot
from mdvtools.charts.histogram_plot import HistogramPlot
from mdvtools.charts.dot_plot import DotPlot
from mdvtools.charts.box_plot import BoxPlot
from mdvtools.charts.scatter_plot_3D import ScatterPlot3D
from mdvtools.charts.row_chart import RowChart
from mdvtools.charts.scatter_plot import ScatterPlot
from mdvtools.charts.abundance_box_plot import AbundanceBoxPlot
from mdvtools.charts.stacked_row_plot import StackedRowChart
from mdvtools.charts.ring_chart import RingChart
from mdvtools.charts.violin_plot import ViolinPlot

ViolinPlot

import json 




def load_data(path):
    #Load data from the specified CSV file.
    return pd.read_csv(path, low_memory=False)

def convert_plot_to_json(plot):
    #Convert plot data to JSON format.
    return json.loads(json.dumps(plot.plot_data, indent=2).replace("\\", ""))
    
def create_heatmap_plot(title, params, size, position, colorscale, x_axis_settings, y_axis_settings):
    """Create and configure a HeatmapPlot instance with the given parameters."""
    plot = HeatmapPlot(
        title=title,
        params=params,
        size=size,
        position=position
    )
    
    plot.set_color_scale(colorscale)  # colorscale setting
    plot.set_x_axis(**x_axis_settings)  # x-axis settings
    plot.set_y_axis(**y_axis_settings)  # y-axis settings
    
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
    data_path = '/Users/petertodd/code/www/MDV-maria/python/sample_data/data_cells.csv'
    view_name = "default"
    
    # Create project
    project = MDVProject(project_path, delete_existing=True)
    
    # Load data
    data_frame = load_data(data_path)

    # Set the correct data type to the "localisation_status" data source (imports as integer but it should be str to appear as a category)
    data_frame['localisation_status'] = data_frame['localisation_status'].apply(str)
    
    # Add datasource
    project.add_datasource(data_path, data_frame)
    
    # HeatmapPlot parameters
    title = "Heatmap Plot of Localisation Status vs UTR Length"
    params = ["localisation_status", "UTR_length"]
    size = [792, 472]
    position = [10, 10]
    
    colorscale = {
        'log': False
    }
    
    x_axis_settings = {
        'axis_labels': "Localisation Status",
        'axis_title': "Localisation Status",
    }
    
    y_axis_settings = {
        'axis_labels': "UTR Length",
        'axis_title': "UTR Length",
    }
    
    # Create plot
    heatmap_plot = create_heatmap_plot(title, params, size, position, colorscale, x_axis_settings, y_axis_settings)

    # Convert plot to JSON and set view
    heatmap_plot_json = convert_plot_to_json(heatmap_plot)
    heatmap_view = {'initialCharts': {data_path: [heatmap_plot_json]}}
    
    project.set_view(view_name, heatmap_view)
    project.set_editable(True)
    project.serve()

if __name__ == "__main__":
    main()