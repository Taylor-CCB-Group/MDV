## This file creates a new MDV project named 'project' containing a view named 'default'
## The datasource used is an h5ad file that was provided at run time
## The view 'default' shows two scatterplots
## The h5ad file is an AnnData object and here the obs attribute was used
## The first scatter plot uses n_genes_by_counts (number of genes detected per cell) versus total_counts (total gene expression counts per cell). 
## This plot provides a quick quality control check to identify cells with low or excessively high gene counts, which may indicate low-quality cells (e.g., empty droplets) or doublets. 
## Such checks help refine the dataset for downstream analysis.
## The second scatter plot displays S_score (S phase) versus G2M_score (G2/M phase), representing cell cycle phases based on gene expression profiles. 
## Cell cycle scores are relevant in single-cell biology as they can impact clustering and differential expression analysis. 
## This visualization helps detect and account for cell cycle effects in the dataset.


import os
import pandas as pd
import scanpy as sc
import sys
import json 
from mdvtools.mdvproject import MDVProject
from mdvtools.charts.scatter_plot import ScatterPlot

    
def create_scatter_plot(title, params, size, position, color, x_axis_settings, y_axis_settings):
    """Create and configure a ScatterPlot instance with the given parameters."""
    plot = ScatterPlot(
        title=title,
        params= params,
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
    data_path = "file_path"
    adata = sc.read_h5ad(data_path)
    cells_df = pd.DataFrame(adata.obs)

    # Name datasource
    datasource_name = 'cells'
    cells_df.name = datasource_name
    
    # Create project
    project = MDVProject(project_path, delete_existing=True)
    
    # Add datasources
    project.add_datasource(datasource_name, cells_df)
    
    # ScatterPlot parameters for the first scatter plot
    title1 = "Scatter Plot 1"
    params = ["n_genes_by_counts", "total_counts"]
    size1 = [792, 472]
    position1 = [10, 10]

    color = 'final_analysis'
    
    x_axis_settings1 = {
        'size': 30,
        'label': "Number of Genes by Counts",
        'textsize': 13,
        'tickfont': 10
    }
    
    y_axis_settings1 = {
        'size': 45,
        'label': "Total Counts",
        'textsize': 13,
        'tickfont': 10,
        'rotate_labels': False
    }
    
    # ScatterPlot parameters for the second scatter plot
    title2 = "Scatter Plot 2"
    params = ["S_score","G2M_score"]
    size2 = [792, 472]
    position2 = [820, 10]
    
    x_axis_settings2 = {
        'size': 30,
        'label': "S Score",
        'textsize': 13,
        'tickfont': 10
    }
    
    y_axis_settings2 = {
        'size': 45,
        'label': "G2M Score",
        'textsize': 13,
        'tickfont': 10,
        'rotate_labels': False
    }
    
    # Create and configure scatter plots
    scatter_plot1 = create_scatter_plot(
        title1, params, size1, position1, color, x_axis_settings1, y_axis_settings1
    )
    
    scatter_plot2 = create_scatter_plot(
        title2, params, size2, position2, color, x_axis_settings2, y_axis_settings2
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