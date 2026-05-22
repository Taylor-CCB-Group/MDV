"""
ChatMDV RAG example (BoxPlot). Use set_axis_properties for axis labels — never set_x_axis / set_y_axis
(those exist only on HistogramPlot, HeatmapPlot, MultiLinePlot, AbundanceBoxPlot).
"""
import os
import json
import pandas as pd
import scanpy as sc
from mdvtools.mdvproject import MDVProject
from mdvtools.charts.box_plot import BoxPlot

def create_box_plot(title, params, size, position, x_axis_settings, y_axis_settings):
    """Create and configure a BoxPlot instance with the given parameters."""
    plot = BoxPlot(
        title=title,
        params=params,
        size=size,
        position=position
    )
    plot.set_axis_properties("x", x_axis_settings)
    plot.set_axis_properties("y", y_axis_settings)
    return plot

def convert_plot_to_json(plot):
    """Convert plot data to JSON format."""
    return json.loads(json.dumps(plot.plot_data, indent=2).replace("\\\\", ""))

def main():
    """Main function to create the project and serve it."""
    # Constants
    project_path = os.path.expanduser('~/mdv/project')
    data_path = "path_to_data"
    view_name = "Example view name"
    datasource_name = "datasource_name"
    datasource_name_2 = "datasource_name_2"
    
    # Create project
    project = MDVProject(project_path, delete_existing=True)
    
    # Load data
    adata = sc.read_h5ad(data_path, backed='r')
    data_frame_obs = pd.DataFrame(adata.obs)

    data_frame_var = pd.DataFrame(adata.var)
    data_frame_var['name'] = adata.var_names.to_list()
    
    # Add datasource
    project.add_datasource(datasource_name, data_frame_obs)
    project.add_datasource(datasource_name_2, data_frame_var)
    
    # BoxPlot parameters
    title = "Example title"
    params = ["param1", "param2"] #param1 should be a categorical variable and param2 should be a numerical variable.
    size = [1092, 472]
    position = [10, 10]
    
    x_axis_settings = {"label": "Category", "textSize": 13, "tickfont": 10}
    y_axis_settings = {"label": "Value", "textSize": 13, "tickfont": 10}

    # Create plot
    plot = create_box_plot(title, params, size, position, x_axis_settings, y_axis_settings)
    
    # Convert plot to JSON and set view
    boxplot_chart_json = convert_plot_to_json(plot)
    boxplot_view = {'initialCharts': {datasource_name: [boxplot_chart_json]}}
    
    project.set_view(view_name, boxplot_view)
    project.set_editable(True)
    project.serve()

if __name__ == "__main__":
    main()
