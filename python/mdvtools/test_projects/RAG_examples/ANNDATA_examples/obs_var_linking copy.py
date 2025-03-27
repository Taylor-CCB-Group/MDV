import os
import pandas as pd
import scanpy as sc
from mdvtools.mdvproject import MDVProject
from mdvtools.charts.dot_plot import DotPlot
import json 

def should_add_link(project, datasource, link_name):
    """
    Returns True if the given link_name or its expected structure is missing or mismatched.
    Returns False if everything is present and correctly structured.
    """
    links = project.get_datasource_metadata(datasource).get("links") 
    if not links:
        return True  # links missing
    link_entry = links.get("genes")
    if not link_entry:
        return True  # 'genes' link missing
    rows_as_columns = link_entry.get("rows_as_columns")
    if not rows_as_columns:
        return True  # rows_as_columns missing
    subgroups = rows_as_columns.get("subgroups")
    if not subgroups:
        return True  # subgroups missing or empty
    name_field = rows_as_columns.get("name")
    if name_field != link_name:
        return True  # name field doesn't match expected link_name
    if link_name not in subgroups:
        return True  # specified link_name not present in subgroups
    return False  # structure is complete and valid



def create_dot_plot(title, params, size, position, colorscale):
    """Create and configure a DotPlot instance with the given parameters."""
    plot = DotPlot(
        title=title,
        params=params,
        size=size,
        position=position
    )
    plot.set_color_scale(colorscale)
    
    return plot

def convert_plot_to_json(plot):
    return json.loads(json.dumps(plot.plot_data, indent=2).replace("\\", ""))

def main():
    """Main function to create the project and serve it."""
        # Constants
    project_path = os.path.expanduser('~/mdv/project')
    data_path = "path_to_data"
    view_name = "default"
    datasource_name = "datasource_name"
    datasource_name_2 = "datasource_name_2"

    # Create project
    project = MDVProject(project_path, delete_existing=True)
    
    # Load data
    adata = sc.read_h5ad(data_path)
    data_frame_obs = pd.DataFrame(adata.obs)

    data_frame_var = pd.DataFrame(adata.var)
    data_frame_var['name'] = adata.var_names.to_list()

    
    # Add datasource
    project.add_datasource(datasource_name, data_frame_obs)
    project.add_datasource(datasource_name_2, data_frame_var)
    
    # Update datasource
    # project.set_column(datasource_name_2, "variable_name", data_frame_var['variable_name'])

    # DotPlot parameters
    param2 = "param2"
    param2_index = data_frame_var['name'].tolist().index(param2)

    # DotPlot parameters
    dot_title = "Dot plot example title"
    dot_params = ["param1", f"link|{param2}(link)|{param2_index}"]
    dot_size = [400, 250]
    dot_position = [10, 500]

    dot_colorscale = {
        'log': False
    }
    
    # Create dot plot
    dot_plot = create_dot_plot(dot_title, dot_params, dot_size, dot_position, dot_colorscale)

    dot_plot_json = convert_plot_to_json(dot_plot)
    
    view_config = {'initialCharts': {datasource_name: [dot_plot_json]}}

    # creating the link between the two datasets so that selecting a subset of genes to add the expression in cells is enabled
    if should_add_link(project, datasource_name, "link"):
        project.add_rows_as_columns_link(datasource_name,datasource_name_2,"name","link")
        project.add_rows_as_columns_subgroup(datasource_name,datasource_name_2,"link",adata.X)

    project.set_view(view_name, view_config)
    project.set_editable(True)
    project.serve()

if __name__ == "__main__":
    main()