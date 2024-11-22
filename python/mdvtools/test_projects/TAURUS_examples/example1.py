## This file creates a new MDV project named 'project' containing a view named 'default'
## The datasource used is an h5ad file that was provided at run time
## The view 'default' shows a histogram of the doublet scores data
## The h5ad file is an AnnData object and here the obs attribute was used
## The doublet score in the obs attribute of an AnnData object is typically a numeric value assigned to each cell 
## to indicate the likelihood that it is a doublet—a single observation that actually represents two (or more) cells 
## that were mistakenly encapsulated together during cell capture.

import os
import pandas as pd
import scanpy as sc
import json 
import sys

from mdvtools.mdvproject import MDVProject
from mdvtools.charts.histogram_plot import HistogramPlot


    
def create_histogram_plot(title, param, bin_number, display_min, display_max, size, position, x_axis_settings, y_axis_settings):
    """Create and configure a HistogramPlot instance with the given parameters."""
    plot = HistogramPlot(
        title=title,
        param=param,  # When the variable is named as "param", it can only take one data field.
        bin_number=bin_number,
        display_min=display_min,
        display_max=display_max,
        size=size,
        position=position
    )
    
    plot.set_x_axis(**x_axis_settings)  # x-axis settings
    plot.set_y_axis(**y_axis_settings)  # y-axis settings
    
    return plot

def convert_plot_to_json(plot):
    """Convert plot data to JSON format."""
    return json.loads(json.dumps(plot.plot_data, indent=2).replace("\\", ""))

def main():
    """Main function to create the project and serve it."""
    # Constants
    project_path = os.path.expanduser('~/mdv/project')
    view_name = "default"
    
    # Load data
    data_path = "file_path"
    adata = sc.read_h5ad(data_path)
    cells_df = pd.DataFrame(adata.obs)

    # Name datasource
    datasource_name = 'cells'
    cells_df.name = datasource_name
    
    # Create project
    project = MDVProject(project_path, delete_existing=True)
    
    # Add datasource
    project.add_datasource(datasource_name, cells_df)
    
    # HistogramPlot parameters
    title = "Doublet Scores Distribution"
    param = "doublet_scores"
    bin_number = 50
    display_min = cells_df[param].min()
    display_max = cells_df[param].max()
    size = [792, 472]
    position = [10, 10]
    
    x_axis_settings = {
        'size': 30,
        'label': "Doublet Scores",
        'textsize': 13,
        'tickfont': 10
    }
    
    y_axis_settings = {
        'size': 45,
        'label': "Frequency",
        'textsize': 13,
        'tickfont': 10,
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