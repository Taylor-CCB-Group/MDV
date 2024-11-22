## This file creates a new MDV project named 'project' containing a view named 'default'
## The datasource used is an h5ad file that was provided at run time
## The view 'default' shows one stacked bar plot.
## The h5ad file is an AnnData object and here the obs attribute was used.
## The stacked bar plot shows the distribution of various cell states across different patients.
## By examining the cell state proportions for each patient, one can observe if there are unique patterns or 
## distributions specific to certain patients. For example, patients with particular conditions may show a higher abundance of certain cell states.


import os
import pandas as pd
import scanpy as sc
import sys
import json
from mdvtools.mdvproject import MDVProject
from mdvtools.charts.stacked_row_plot import StackedRowChart
    
def create_stacked_row_plot(title, params, size, position, legend_display, xaxis_properties, yaxis_properties):
    """Create and configure a StackedRowChart instance with the given parameters."""
    plot = StackedRowChart(
        title=title,
        params=params,
        size=size,
        position=position
    )

    plot.set_color_legend(legend_display)
    plot.set_axis_properties("x", xaxis_properties)
    plot.set_axis_properties("y", yaxis_properties)

    return plot

def convert_plot_to_json(plot):
    """Convert plot data to JSON format."""
    return json.loads(json.dumps(plot.plot_data, indent=2).replace("\\\\", ""))

def main():
    """Main function to create the project and serve it."""
    # Constants
    project_path = os.path.expanduser('~/mdv/project')
    view_name = "default"
    
    # Load data
    data_path = "file_path"
    adata = sc.read_h5ad(data_path)
    cells_df = pd.DataFrame(adata.obs)
    cells_df.name = 'cells'
    
    # Create project
    project = MDVProject(project_path, delete_existing=True)
    
    # Add datasource
    project.add_datasource('cells', cells_df)
    
    # StackedRowChart parameters
    title = "Abundance of Cell Types per Patient"
    params = ["final_analysis", "Patient"]  # Using 'final_analysis' for cell state and 'Patient' for patient
    size = [792, 472]
    position = [10, 10]

    legend_display = True
    
    xaxis_properties = {"label": "Patient", 
                        "textSize": 13, 
                        "tickfont": 10
    }

    yaxis_properties = {"label": "Cell Type", 
                        "textSize": 13, 
                        "tickfont": 10
    }

    # Create plot
    stacked_row_plot = create_stacked_row_plot(title, params, size, position, legend_display, xaxis_properties, yaxis_properties)
    
    # Convert plot to JSON and set view
    stacked_row_plot_json = convert_plot_to_json(stacked_row_plot)
    stackedrowchart_view = {'initialCharts': {'cells': [stacked_row_plot_json]}}
    
    project.set_view(view_name, stackedrowchart_view)
    project.set_editable(True)
    project.serve()

if __name__ == "__main__":
    main()