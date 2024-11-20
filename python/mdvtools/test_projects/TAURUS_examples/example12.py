## This file creates a new MDV project named 'project' containing a view named 'default'
## The datasource used is an h5ad file that was provided at run time
## The view 'default' shows one box plot.
## The h5ad file is an AnnData object and here the obs attribute was used.
## The format f"Gene expression|{gene_name}(Gene expression)|{genes_df.index.get_loc(gene_name)}" is used to get the gene expression data for that specific gene
## By visualizing MIR1302-2HG expression across diseases, researchers can explore how its expression levels differ between disease states and potentially identify patterns. 
## For instance, consistently higher MIR1302-2HG levels in specific conditions may confirm its role as a disease marker or therapeutic target.



import os
import pandas as pd
import scanpy as sc
import sys
from mdvtools.mdvproject import MDVProject
from mdvtools.charts.box_plot import BoxPlot
import json 


def create_box_plot(title, params, size, position):
    """Create and configure a BoxPlot instance with the given parameters."""
    plot = BoxPlot(
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
    adata = sc.read_h5ad(sys.argv[1])
    cells_df = pd.DataFrame(adata.obs)
    cells_df.name = 'cells'
    
    genes_df = pd.DataFrame(adata.var)
    genes_df.name = 'genes'
    genes_df['gene_id'] = genes_df.index

    # Rename 'final_analysis' to 'cell state'
    cells_df.rename(columns={"final_analysis": "cell state"}, inplace=True)
    
    # Create project
    project = MDVProject(project_path, delete_existing=True)
    
    # Add datasource
    project.add_datasource('cells', cells_df)
    project.add_datasource('genes', genes_df)
    
    # Create a link between the two datasets
    project.add_rows_as_columns_link("cells", "genes", "gene_id", "Gene Expression")
    project.add_rows_as_columns_subgroup("cells", "genes", "Gene expression", adata.X.toarray())
    
    # BoxPlot parameters for the specific gene "MIR1302-2HG"
    gene_name = "MIR1302-2HG"

    # The format f"Gene expression|{gene_name}(Gene expression)|{genes_df.index.get_loc(gene_name)}" is used to get the gene expression data for that specific gene
    box_title = f"Gene expression for {gene_name} per Disease"
    box_params = ["Disease", f"Gene expression|{gene_name}(Gene expression)|{genes_df.index.get_loc(gene_name)}"]
    box_size = [792, 472]
    box_position = [10, 10]
    
    # Create plot
    plot = create_box_plot(box_title, box_params, box_size, box_position)
    
    # Convert plot to JSON and set view
    boxplot_chart_json = convert_plot_to_json(plot)
    boxplot_view = {'initialCharts': {'cells': [boxplot_chart_json]}}
    
    project.set_view(view_name, boxplot_view)
    project.set_editable(True)
    project.serve()

if __name__ == "__main__":
    main()