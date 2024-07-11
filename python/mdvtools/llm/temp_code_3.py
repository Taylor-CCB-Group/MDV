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
    
def create_scatter_plot(title, params, size, position, color, brush, opacity, radius, legend_display, legend_position, xaxis_properties, yaxis_properties):
    """Create and configure a ScatterPlot instance with the given parameters."""
    plot = ScatterPlot(
        title=title,
        params=params,
        size=size,
        position=position
    )

    plot.set_default_color(color)
    plot.set_brush(brush)
    plot.set_opacity(opacity)
    plot.set_radius(radius)
    plot.set_color_legend(legend_display, legend_position)
    plot.set_axis_properties("x", xaxis_properties)
    plot.set_axis_properties("y", yaxis_properties)
    
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
    data_path = '/Users/petertodd/code/www/MDV-maria/python/mdvtools/llm/sample_data/data_cells.csv'
    view_name = "default"
    
    # Create project
    project = MDVProject(project_path, delete_existing=True)
    
    # Load data
    data_frame = load_data(data_path)

    # Set the correct data type to the "leiden" data source (imports as integer but it should be str to appear as a category)
    data_frame['leiden'] = data_frame['leiden'].apply(str)
    
    # Add datasource
    project.add_datasource(data_path, data_frame)
    
    # ScatterPlot parameters
    title = "Scatter Plot of UMAP Coloured by Cluster"
    params = ["X_umap_1", "X_umap_2", "leiden"]
    size = [792, 472]
    position = [10, 10]

    color = "#377eb8"
    brush = "poly"
    opacity = 0.8
    radius = 0.2

    legend_display = True
    legend_position = [375,1]
              
    xaxis_properties = {"label": "X_umap_1", 
             "textSize": 13, 
             "tickfont": 10
    }

    yaxis_properties = {"label": "X_umap_2", 
             "textSize": 13, 
             "tickfont": 10
    }

    # Create plot
    scatter_plot = create_scatter_plot(title, params, size, position, color, brush, opacity, radius, legend_display, legend_position, xaxis_properties, yaxis_properties)
    
    # Convert plot to JSON and set view
    scatter_plot_json = convert_plot_to_json(scatter_plot)
    scatterplot_view = {'initialCharts': {data_path: [scatter_plot_json]}}
    
    project.set_view(view_name, scatterplot_view)
    project.set_editable(True)
    project.serve()

if __name__ == "__main__":
    main()