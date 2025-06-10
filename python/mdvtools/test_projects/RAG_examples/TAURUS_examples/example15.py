import os
import pandas as pd
import scanpy as sc
from mdvtools.mdvproject import MDVProject
from mdvtools.charts.dot_plot import DotPlot
import json
import numpy as np    


def create_dot_plot(title, params, size, position):
    """Create and configure a DotPlot instance with the given parameters."""
    plot = DotPlot(
        title=title,
        params=params,
        size=size,
        position=position
    )
    plot.set_axis_properties("x", {"label": "Cell State", "textSize": 13, "tickfont": 10})
    plot.set_axis_properties("y", {"label": "Expression", "textSize": 13, "tickfont": 10})
    plot.set_color_scale(log_scale=False)
    plot.set_color_legend(True, [40, 10])
    plot.set_fraction_legend(True, [0, 0])
    return plot

def convert_plot_to_json(plot):
    """Convert plot data to JSON format."""
    return json.loads(json.dumps(plot.plot_data, indent=2).replace("\\", ""))

def main():
    """Main function to create the project and serve it."""
    # Constants
    project_path = os.path.expanduser('~/mdv/project')
    view_name = "Gene with the lowest expression, can I see it's dot plot per cell state?"
    
    # Load data
    data_path = "file_path"
    adata = sc.read_h5ad(data_path)
    cells_df = pd.DataFrame(adata.obs)
    genes_df = pd.DataFrame(adata.var)
    genes_df['gene_id'] = genes_df.index

    # Find the gene with the lowest expression
    gene_expression_means = adata.X.mean(axis=0) # type: ignore
    lowest_expression_index = np.argmin(gene_expression_means)
    lowest_expression_gene = genes_df.index[lowest_expression_index] # type: ignore

   # Create project
    project = MDVProject(project_path, delete_existing=True)
    
    # Add datasource
    project.add_datasource('cells', cells_df)
    project.add_datasource('genes', genes_df)
    
    # Create a link between the two datasets
    project.add_rows_as_columns_link("cells", "genes", "gene_id", "Gene Expression")
    project.add_rows_as_columns_subgroup("cells", "genes", "Gene expression", adata.X)
    
    # DotPlot parameters for the gene with the lowest expression
    dot_title = f"Gene expression for {lowest_expression_gene} per cell state"
    dot_params = ['sample_id',  'final_analysis', f'Gene expression|{lowest_expression_gene}(Gene expression)|{lowest_expression_index}']
    dot_size = [450, 300]
    dot_position = [10, 10]
    
    # Create dot plot
    dot_plot = create_dot_plot(dot_title, dot_params, dot_size, dot_position)
    
    # Convert plot to JSON and set view
    dot_plot_json = convert_plot_to_json(dot_plot)
    view_config = {'initialCharts': {'cells': [dot_plot_json]}}
    
    project.set_view(view_name, view_config)
    project.set_editable(True)
    # project.serve()

if __name__ == "__main__":
    main()
else:
    main()