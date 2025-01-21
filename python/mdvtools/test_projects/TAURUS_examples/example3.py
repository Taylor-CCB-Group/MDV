## This file creates a new MDV project named 'project' containing a view named 'default'
## The datasource used is an h5ad file that was provided at run time
## The view 'default' shows two dotplots
## The h5ad file is an AnnData object and here the obs attribute was used
## The first dot plot visualizes total counts (often corresponding to total gene expression counts per cell) across different patients.
## It helps identify inter-patient variability in total expression levels.
## For example, variation here could reflect differences in sample preparation, patient biology, or experimental noise.
## The second dot plot displays total counts across different disease states.
## This could highlight how overall cell expression levels vary by disease status, which might indicate disease-related cellular activation or suppression patterns.

import os
import pandas as pd
import scanpy as sc
import sys
import json
from mdvtools.mdvproject import MDVProject
from mdvtools.charts.dot_plot import DotPlot

    
def create_dot_plot(title, params, size, position, x_axis_settings, y_axis_settings):
    """Create and configure a DotPlot instance with the given parameters."""
    plot = DotPlot(
        title=title,
        params=params,
        size=size,
        position=position
    )

    plot.set_axis_properties("x", x_axis_settings)  # x-axis settings
    plot.set_axis_properties("y", y_axis_settings)  # y-axis settings

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
    datasource_name = "datasource_name"
    
    # Create project
    project = MDVProject(project_path, delete_existing=True)
    
    # Add datasource
    project.add_datasource(datasource_name, cells_df)
    
    # DotPlot parameters for the first dot plot
    title1 = "Dot Plot 1"
    params1 = ["Patient", "total_counts"]
    size1 = [792, 472]
    position1 = [10, 10]
    
    x_axis_settings1 = {
        'size': 30,
        'label': "Total Counts",
        'textsize': 13,
        'tickfont': 10
    }
    
    y_axis_settings1 = {
        'size': 45,
        'label': "Patient", 
        'textsize': 13,
        'tickfont': 10,
        'rotate_labels': False
    }
    
    # DotPlot parameters for the second dot plot
    title2 = "Dot Plot 2"
    params2 = ["Disease", "total_counts"]
    size2 = [792, 472]
    position2 = [820, 10]
    
    x_axis_settings2 = {
        'size': 30,
        'label': "Total counts",
        'textsize': 13,
        'tickfont': 10
    }
    
    y_axis_settings2 = {
        'size': 45,
        'label': "Disease",
        'textsize': 13,
        'tickfont': 10,
        'rotate_labels': False
    }
    
    # Create and configure dot plots
    dot_plot1 = create_dot_plot(
        title1, params1, size1, position1, x_axis_settings1, y_axis_settings1
    )
    
    dot_plot2 = create_dot_plot(
        title2, params2, size2, position2, x_axis_settings2, y_axis_settings2
    )
    
    # Convert plots to JSON and set view
    dot_chart_json1 = convert_plot_to_json(dot_plot1)
    dot_chart_json2 = convert_plot_to_json(dot_plot2)
    
    dot_view = {'initialCharts': {datasource_name: [dot_chart_json1, dot_chart_json2]}}
    
    project.set_view(view_name, dot_view)
    project.set_editable(True)
    project.serve()

if __name__ == "__main__":
    main()