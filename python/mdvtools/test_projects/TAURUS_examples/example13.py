import os
import sys
import pandas as pd
import scanpy as sc
import  numpy as np
from mdvtools.mdvproject import MDVProject
from mdvtools.charts.dot_plot import DotPlot
import json 


def create_dot_plot(title, params, size, position):
    plot = DotPlot(
        title=title,
        params=params,
        size=size,
        position=position
    )
    plot.set_axis_properties("x", {"label": "", "textSize": 13, "tickfont": 10})
    plot.set_axis_properties("y", {"label": "", "textSize": 13, "tickfont": 10})
    plot.set_color_scale(log_scale=False)
    plot.set_color_legend(True, [40, 10])
    plot.set_fraction_legend(True, [0, 0])
    return plot

def convert_plot_to_json(plot):
    return json.loads(json.dumps(plot.plot_data, indent=2).replace("\\", ""))

def main():
    project_path = os.path.expanduser('~/mdv/project')
    view_name = "default"
    
    # Load data
    data_path = "file_path"
    adata = sc.read_h5ad(data_path)
    cells_df = pd.DataFrame(adata.obs)
    cells_df.name = 'cells'
    
    genes_df = pd.DataFrame(adata.var)
    genes_df['gene_id'] = genes_df.index
    
    # Rename 'final_analysis' to 'cell state'
    cells_df.rename(columns={"final_analysis": "cell state"}, inplace=True)
    
    # Add UMAP data to the dataframe
    umap_np = np.array(adata.obsm["X_umap"])
    cells_df["UMAP 1"] = umap_np[:, 0]
    cells_df["UMAP 2"] = umap_np[:, 1]
    
    # Create project
    project = MDVProject(project_path, delete_existing=True)
    
    # Add datasource
    project.add_datasource('cells', cells_df)
    project.add_datasource('genes', genes_df)
    
    # DotPlot parameters for CD13 (ANPEP)
    gene_name = "ANPEP"
    gene_id = "ENSG00000166825"
    dot_title = f"Gene expression for {gene_name} per cell state"
    dot_params = ["cell state", f"Gene expression|{gene_name}(Gene expression)|{genes_df.index.get_loc(gene_name)}"]
    dot_size = [450, 300]
    dot_position = [10, 10]
    
    # Create dot plot
    dot_plot = create_dot_plot(dot_title, dot_params, dot_size, dot_position)
    
    # Convert plots to JSON and set view
    dot_plot_json = convert_plot_to_json(dot_plot)
    
    view_config = {'initialCharts': {'cells': [dot_plot_json]}}
    
    # Creating the link between the two datasets so that selecting a subset of genes to add the expression in cells is enabled
    project.add_rows_as_columns_link("cells", "genes", "gene_id", "Gene expression")
    project.add_rows_as_columns_subgroup("cells", "genes", "Gene expression", adata.X.toarray()) # add the gene expression
    
    project.set_view(view_name, view_config)
    project.set_editable(True)
    project.serve()

if __name__ == "__main__":
    main()