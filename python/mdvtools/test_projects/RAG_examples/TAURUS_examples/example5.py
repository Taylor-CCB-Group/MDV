## This file creates a new MDV project named 'project' containing a view named 'default'
## The datasource used is an h5ad file that was provided at run time
## The view 'default' shows one table plot
## The h5ad file is an AnnData object and here the obs attribute was used
## The table allows viewing of all important cell features in one place, making it easy to sort, or review cell-specific metrics across 
## different categories like patient ID, disease state, or treatment.

import os
import pandas as pd
import scanpy as sc
from mdvtools.mdvproject import MDVProject
from mdvtools.charts.table_plot import TablePlot
import json 




def load_data(path):
    #Load data from the specified CSV file.
    return pd.read_csv(path, low_memory=False)

def convert_plot_to_json(plot):
    #Convert plot data to JSON format.
    return json.loads(json.dumps(plot.plot_data, indent=2).replace("\\", ""))
    
def create_table_plot(title, params, size, position):
    """Create and configure a TablePlot instance with the given parameters."""
    plot = TablePlot(
        title=title,
        params=params,
        size=size,
        position=position
    )
    
    return plot

def convert_plot_to_json(plot):
    """Convert plot data to JSON format."""
    return json.loads(json.dumps(plot.plot_data, indent=2).replace("\\\\", ""))

def main():
    """Main function to create the project and serve it."""
    # Constants
    project_path = os.path.expanduser('~/mdv/project')
    view_name = "default"
    
    # Load data using scanpy
    data_path = "file_path"
    adata = sc.read_h5ad(data_path)
    cells_df = pd.DataFrame(adata.obs)
    datasource_name = "datasource_name"
    
    # Create project
    project = MDVProject(project_path, delete_existing=True)
    
    # Add datasource
    project.add_datasource(datasource_name, cells_df)

    # TablePlot parameters
    title = "All Data Table"
    params = ['sample_id', 'doublet_scores', 'predicted_doublets', 'n_genes_by_counts', 'total_counts', 
              'total_counts_mt', 'pct_counts_mt', 'total_counts_rp', 'pct_counts_rp', 'total_counts_hb', 
              'pct_counts_hb', 'total_counts_ig', 'pct_counts_ig', 'S_score', 'G2M_score', 'phase', 
              'Disease', 'Patient', 'Site', 'sub_bucket', 'final_analysis', 'Treatment', 'MM_scaled', 
              'Inflammation', 'Remission']
    size = [792, 472]
    position = [10, 10]

    # Create plot
    table_plot = create_table_plot(title, params, size, position)
    
    # Convert plot to JSON and set view
    table_plot_json = convert_plot_to_json(table_plot)
    tableplot_view = {'initialCharts': {datasource_name: [table_plot_json]}}
    
    project.set_view(view_name, tableplot_view)
    project.set_editable(True)
    project.serve()

if __name__ == "__main__":
    main()