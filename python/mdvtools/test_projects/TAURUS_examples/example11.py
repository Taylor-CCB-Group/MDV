## This file creates a new MDV project named 'project' containing a view named 'default'
## The datasource used is an h5ad file that was provided at run time
## The view 'default' shows one box plot.
## The h5ad file is an AnnData object and here the obs attribute was used.
## The format f"Gene expression|{gene_name}(Gene expression)|{genes_df.index.get_loc(gene_name)}" is used to get the gene expression data for that specific gene
## By visualizing TNF expression across diseases, researchers can explore how its expression levels differ between disease states and potentially identify patterns. 
## For instance, consistently higher TNF levels in specific conditions may confirm its role as a disease marker or therapeutic target.

import os
import pandas as pd
import scanpy as sc
import sys
from mdvtools.mdvproject import MDVProject
from mdvtools.charts.box_plot import BoxPlot
import json 
    
def create_box_plot(title, params, size, position, plot_id):
    """Create and configure a BoxPlot instance with the given parameters."""
    plot = BoxPlot(
        title=title,
        params=params,
        size=size,
        position=position,
        id=plot_id
    )
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
    
    genes_df = pd.DataFrame(adata.var)
    genes_df['gene_id'] = genes_df.index
    
    # Create project
    project = MDVProject(project_path, delete_existing=True)
    
    # Add datasource
    project.add_datasource('cells', cells_df)
    project.add_datasource('genes', genes_df)
    
    # BoxPlot parameters
    gene_name = "TNF"  # Choose a gene with interesting properties

    # The format f"Gene expression|{gene_name}(Gene expression)|{genes_df.index.get_loc(gene_name)}" is used to get the gene expression data for that specific gene
    box_title = f"Gene expression for {gene_name} per disease box plot"
    box_params = ["Disease", f"Gene expression|{gene_name}(Gene expression)|10452"]
    box_size = [615, 557]
    box_position = [50, 50]
    box_plot_id = "boxPlot1"
    
    # Create box plot
    box_plot = create_box_plot(
        box_title, box_params, box_size, box_position, box_plot_id
    )
    
    # Convert plot to JSON and set view
    box_plot_json = convert_plot_to_json(box_plot)
    view_config = {'initialCharts': {'cells': [box_plot_json]}}
    
    # Create the link between the two datasets
    project.add_rows_as_columns_link("cells", "genes", "gene_id", "Gene Expression")
    project.add_rows_as_columns_subgroup("cells", "genes", "Gene expression", adata.X.toarray())  # Add the gene expression
    
    project.set_view(view_name, view_config)
    project.set_editable(True)
    project.serve()

if __name__ == "__main__":
    main()