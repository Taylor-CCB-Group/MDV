import os
import pandas as pd
import scanpy as sc
from mdvtools.mdvproject import MDVProject
from mdvtools.charts.heatmap_plot import HeatmapPlot
import json 

def create_heatmap(title, params, size, position, xaxis_properties, yaxis_properties):
    plot = HeatmapPlot(
        title=title,
        params=params,
        size=size,
        position=position
    )
    plot.set_axis_properties("x", xaxis_properties)
    plot.set_axis_properties("y", yaxis_properties)
    return plot

def convert_plot_to_json(plot):
    return json.loads(json.dumps(plot.plot_data, indent=2).replace("\\", ""))

def main():
    project_path = os.path.expanduser('~/mdv/project')
    view_name = "Heatmap Plot Example"
    
    # Load data
    data_path = "file_path"
    adata = sc.read_h5ad(data_path)
    cells_df = pd.DataFrame(adata.obs)
    datasource_name = "datasource_name"
    
    # Create project
    project = MDVProject(project_path, delete_existing=True)
    
    # Add datasource
    project.add_datasource(datasource_name, cells_df)
    
    # Heatmap parameters
    heatmap_title = "Heatmap: Gene Expression"
    heatmap_params = ['final_analysis','Disease_duration', 'Age', 'Inflammation_score']
    heatmap_size = [792, 472]
    heatmap_position = [10, 10]
    heatmap_xaxis_properties = {"label": "Cells", "textSize": 13, "tickfont": 10}
    heatmap_yaxis_properties = {"label": "Features", "textSize": 13, "tickfont": 10}
    
    # Create heatmap plot
    heatmap_plot = create_heatmap(
        heatmap_title, heatmap_params, heatmap_size, heatmap_position, heatmap_xaxis_properties, heatmap_yaxis_properties
    )
    
    # Convert plot to JSON and set view
    heatmap_plot_json = convert_plot_to_json(heatmap_plot)
    
    view_config = {'initialCharts': {datasource_name: [heatmap_plot_json]}}
    
    project.set_view(view_name, view_config)
    project.set_editable(True)
    project.serve()

if __name__ == "__main__":
    main()