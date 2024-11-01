## This file creates a new MDV project named 'project' containing a view named 'default'
## The datasource used is an h5ad file that was provided at run time
## The view 'default' shows one dotplot and one scatterplot
## The h5ad file is an AnnData object and here the obs attribute was used
## The dot plot helps in assessing whether there is a relationship between disease status and doublet scores. 
## High doublet scores in certain disease groups might indicate potential issues with sample preparation or distinct cellular interactions related to the disease.
## Understanding the distribution of doublet scores across diseases is essential for quality control and can hint at whether particular disease states are associated with higher or lower rates of cellular multiplets.
## The scatter plot visualizes the relationship between gene diversity and total expression levels, helping to identify different cellular states or types across disease conditions.
## For example, disease states associated with higher or lower total counts or gene diversity might suggest differences in cellular activity, metabolic states, or cell types. 
## Identifying clusters or trends could help in discovering disease-specific cell profiles or shifts in cellular populations.


import os
import pandas as pd
import scanpy as sc
import sys
from mdvtools.mdvproject import MDVProject
from mdvtools.charts.dot_plot import DotPlot
from mdvtools.charts.scatter_plot import ScatterPlot
import json 

    
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

def create_scatter_plot(title, params, size, position, color, x_axis_settings, y_axis_settings):
    """Create and configure a ScatterPlot instance with the given parameters."""
    plot = ScatterPlot(
        title=title,
        params=params,
        size=size,
        position=position
    )

    plot.set_color_by(color)
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
    adata = sc.read_h5ad(sys.argv[1])
    cells_df = pd.DataFrame(adata.obs)
    cells_df.name = 'cells'
    
    # Create project
    project = MDVProject(project_path, delete_existing=True)
    
    # Add datasource
    project.add_datasource('cells', cells_df)
    
    # DotPlot parameters
    dot_title = "Dot Plot Example"
    dot_params = ["Disease", "doublet_scores"]
    dot_size = [792, 472]
    dot_position = [10, 10]
    
    dot_x_axis_settings = {
        'size': 30,
        'label': "Doublet Scores",
        'textsize': 13,
        'tickfont': 10
    }
    
    dot_y_axis_settings = {
        'size': 45,
        'label': "Disease",
        'textsize': 13,
        'tickfont': 10,
        'rotate_labels': False
    }
    
    # ScatterPlot parameters
    scatter_title = "Scatter Plot Example"
    scatter_params = ["n_genes_by_counts", "total_counts"]
    scatter_size = [792, 472]
    scatter_position = [820, 10]
    
    scatter_color = 'Disease'
    
    scatter_x_axis_settings = {
        'size': 30,
        'label': "Number of Genes by Counts",
        'textsize': 13,
        'tickfont': 10
    }
    
    scatter_y_axis_settings = {
        'size': 45,
        'label': "Total Counts",
        'textsize': 13,
        'tickfont': 10,
        'rotate_labels': False
    }
    
    # Create and configure plots
    dot_plot = create_dot_plot(
        dot_title, dot_params, dot_size, dot_position, dot_x_axis_settings, dot_y_axis_settings
    )
    
    scatter_plot = create_scatter_plot(
        scatter_title, scatter_params, scatter_size, scatter_position, scatter_color, scatter_x_axis_settings, scatter_y_axis_settings
    )
    
    # Convert plots to JSON and set view
    dot_plot_json = convert_plot_to_json(dot_plot)
    scatter_plot_json = convert_plot_to_json(scatter_plot)
    
    view_config = {'initialCharts': {'cells': [dot_plot_json, scatter_plot_json]}}
    
    project.set_view(view_name, view_config)
    project.set_editable(True)
    project.serve()

if __name__ == "__main__":
    main()