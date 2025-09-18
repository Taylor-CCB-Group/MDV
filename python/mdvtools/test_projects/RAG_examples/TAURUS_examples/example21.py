import os
import pandas as pd
import scanpy as sc
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
    plot.set_axis_properties("x", {"label": "Cell State", "textSize": 13, "tickfont": 10})
    plot.set_axis_properties("y", {"label": "Percentage of Cells Expressing TNF", "textSize": 13, "tickfont": 10})
    plot.set_color_scale(log_scale=False)
    plot.set_color_legend(True, [40, 10])
    plot.set_fraction_legend(True, [0, 0])
    return plot

def convert_plot_to_json(plot):
    return json.loads(json.dumps(plot.plot_data, indent=2).replace("\\", ""))

def main():
    project_path = os.path.expanduser('~/mdv/project')
    view_name = "TNF Gene Expression Percentage"
    
    # Load data
    data_path = "file_path"
    adata = sc.read_h5ad(data_path)
    cells_df = pd.DataFrame(adata.obs)
    cells_df.name = 'cells'
    
    genes_df = pd.DataFrame(adata.var)
    genes_df['gene_id'] = genes_df.index
    
    # Calculate the percentage of cells expressing TNF
    gene_name = "TNF"
    gene_index = genes_df.index.get_loc(gene_name)
    gene_expression_matrix = adata.X
    tnf_expression = gene_expression_matrix[:, gene_index] # type: ignore
    cells_expressing_tnf = (tnf_expression > 0).sum() # type: ignore
    total_cells = len(cells_df)
    percentage_expressing_tnf = (cells_expressing_tnf / total_cells) * 100
    print(f"Percentage of cells expressing TNF: {percentage_expressing_tnf:.2f}%")
    
    # Create project
    project = MDVProject(project_path, delete_existing=True)
    
    # Add datasource
    project.add_datasource('cells', cells_df)
    project.add_datasource('genes', genes_df)
    
    # Create a link between the two datasets
    project.add_rows_as_columns_link("cells", "genes", "gene_id", "Gene Expression")
    project.add_rows_as_columns_subgroup("cells", "genes", "Gene expression", gene_expression_matrix)
    
    # DotPlot parameters for TNF
    dot_title = f"Percentage of Cells Expressing {gene_name}"
    dot_params = ["final_analysis", f"Gene expression|{gene_name}(Gene expression)|{gene_index}"]
    dot_size = [450, 300]
    dot_position = [10, 10]
    
    # Create dot plot
    dot_plot = create_dot_plot(dot_title, dot_params, dot_size, dot_position)
    
    # Convert plot to JSON and set view
    dot_plot_json = convert_plot_to_json(dot_plot)
    view_config = {'initialCharts': {'cells': [dot_plot_json]}}
    
    project.set_view(view_name, view_config)
    project.set_editable(True)
    project.serve()

if __name__ == "__main__":
    main()