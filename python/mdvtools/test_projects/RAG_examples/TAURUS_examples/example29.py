import os
import pandas as pd
import scanpy as sc
from mdvtools.mdvproject import MDVProject
from mdvtools.charts.stacked_row_plot import StackedRowChart
import json 

def create_stacked_row_plot(title, params, size, position, color_legend, xaxis_properties, yaxis_properties):
    plot = StackedRowChart(
        title=title,
        params=params,
        size=size,
        position=position
    )
    plot.set_color_legend(color_legend["display"], color_legend["pos"])
    plot.set_axis_properties("x", xaxis_properties)
    plot.set_axis_properties("y", yaxis_properties)
    return plot

def convert_plot_to_json(plot):
    return json.loads(json.dumps(plot.plot_data, indent=2).replace("\\", ""))

def main():
    project_path = os.path.expanduser('~/mdv/project')
    view_name = "Cell Distribution in Healthy vs Diseased Samples"
    
    # Load data
    data_path = "file_path"
    adata = sc.read_h5ad(data_path)
    cells_df = pd.DataFrame(adata.obs)
    datasource_name = "datasource_name"
    
    # Create project
    project = MDVProject(project_path, delete_existing=True)
    
    # Add datasource
    project.add_datasource(datasource_name, cells_df)
    
    # StackedRowChart parameters
    stacked_title = "Number of Cells in Healthy vs Diseased Samples"
    stacked_params = ["Disease", "final_analysis"]
    stacked_size = [792, 472]
    stacked_position = [10, 10]
    stacked_color_legend = {"display" : True,
                    "pos" : [375,1]}
    stacked_xaxis_properties = {"label": "Disease", "textSize": 13, "tickfont": 10}
    stacked_yaxis_properties = {"label": "Cell State", "textSize": 13, "tickfont": 10}
    
    # Create stacked row plot
    stacked_row_plot = create_stacked_row_plot(
        stacked_title, stacked_params, stacked_size, stacked_position, stacked_color_legend, stacked_xaxis_properties, stacked_yaxis_properties
    )
    
    # Convert plot to JSON and set view
    stacked_row_plot_json = convert_plot_to_json(stacked_row_plot)
    
    view_config = {'initialCharts': {datasource_name: [stacked_row_plot_json]}}
    
    project.set_view(view_name, view_config)
    project.set_editable(True)
    project.serve()

if __name__ == "__main__":
    main()