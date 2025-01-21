import os
import pandas as pd
import scanpy as sc
from mdvtools.mdvproject import MDVProject
from mdvtools.charts.table_plot import TablePlot
import json 

def convert_plot_to_json(plot):
    """Convert plot data to JSON format."""
    return json.loads(json.dumps(plot.plot_data, indent=2).replace("\\\\", ""))

def create_table_plot(title, params, size, position):
    """Create and configure a TablePlot instance with the given parameters."""
    plot = TablePlot(
        title=title,
        params=params,
        size=size,
        position=position
    )
    return plot

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
    
    # Create a summary table of cell counts per cell type
    cell_type_counts = cells_df['final_analysis'].value_counts().reset_index()
    cell_type_counts.columns = ['Cell Type', 'Number of Cells']
    
    # Create project
    project = MDVProject(project_path, delete_existing=True)
    
    # Add datasource
    project.add_datasource(datasource_name, cell_type_counts)

    # TablePlot parameters
    title = "Cell Type Counts"
    params = ['Cell Type', 'Number of Cells']
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