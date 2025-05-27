## This file creates a new MDV project named 'project' containing a view named 'default'
## The datasource used is an h5ad file that was provided at run time
## The view 'default' shows two scatterplots
## The h5ad file is an AnnData object and here the obs attribute was used
## The first scatter plot uses total_counts (total RNA counts per cell) versus pct_counts_mt (percentage of mitochondrial RNA). 
## This plot distinguishes between healthy cells (moderate 'total_counts' and low 'pct_counts_mt') and low-quality cells (low 'total_counts', high 'pct_counts_mt')
## Such checks help in quality control by identifying potential outliers or problematic cells.
## The second scatter plot displays n_genes_by_counts (number of genes detected per cell) versus total_counts (total RNA counts per cell).
## Points with high 'total_counts' but low 'n_genes_by_counts' may indicate doublets or outliers.
## This visualization helps assess the diversity of gene expression in relation to RNA abundance, providing clues about cell type complexity or technical artifacts.


import os
import pandas as pd
import scanpy as sc
import sys
import json 
from mdvtools.mdvproject import MDVProject
from mdvtools.charts.scatter_plot import ScatterPlot

    
def create_scatter_plot(title, params, size, position, x_axis_settings, y_axis_settings):
    """Create and configure a ScatterPlot instance with the given parameters."""
    plot = ScatterPlot(
        title=title,
        params= params,
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
    
    # Add datasources
    project.add_datasource(datasource_name, cells_df)
    
    # ScatterPlot parameters for the first scatter plot
    title1 = "Scatter Plot 1"
    params = ["total_counts", "pct_counts_mt"]
    size1 = [792, 472]
    position1 = [10, 10]
    
    x_axis_settings1 = {
        'size': 30,
        'label': "Total RNA counts per cell",
        'textsize': 13,
        'tickfont': 10
    }
    
    y_axis_settings1 = {
        'size': 45,
        'label': "Percentage of mitochondrial RNA",
        'textsize': 13,
        'tickfont': 10,
        'rotate_labels': False
    }
    
    # ScatterPlot parameters for the second scatter plot
    title2 = "Scatter Plot 2"
    params = ["n_genes_by_counts","total_counts"]
    size2 = [792, 472]
    position2 = [820, 10]
    
    x_axis_settings2 = {
        'size': 30,
        'label': "Number of genes detected per cell",
        'textsize': 13,
        'tickfont': 10
    }
    
    y_axis_settings2 = {
        'size': 45,
        'label': "Total RNA counts per cell",
        'textsize': 13,
        'tickfont': 10,
        'rotate_labels': False
    }
    
    # Create and configure scatter plots
    scatter_plot1 = create_scatter_plot(
        title1, params, size1, position1, x_axis_settings1, y_axis_settings1
    )
    
    scatter_plot2 = create_scatter_plot(
        title2, params, size2, position2, x_axis_settings2, y_axis_settings2
    )
    
    # Convert plots to JSON and set view
    scatter_chart_json1 = convert_plot_to_json(scatter_plot1)
    scatter_chart_json2 = convert_plot_to_json(scatter_plot2)
    
    scatter_view = {'initialCharts': {datasource_name: [scatter_chart_json1, scatter_chart_json2]}}
    
    project.set_view(view_name, scatter_view)
    project.set_editable(True)
    project.serve()

if __name__ == "__main__":
    main()