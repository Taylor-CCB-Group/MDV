import os
import pandas as pd
import scanpy as sc
from mdvtools.mdvproject import MDVProject
from mdvtools.charts.stacked_row_plot import StackedRowChart
from mdvtools.charts.selection_dialog_plot import SelectionDialogPlot

import json

def create_stacked_row_plot(title, params, size, position, color_legend, xaxis_properties, yaxis_properties):
    """Create and configure a StackedRowChart instance with the given parameters."""
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

def create_selection_dialog_plot(title, params, size, position):
    """Create and configure a SelectionDialogPlot instance with the given parameters."""
    plot = SelectionDialogPlot(
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
    
    # StackedRowChart parameters
    stacked_row_title = "Stacked Row Plot Example"
    stacked_row_params = ["param1", "param3"]  # Use two categorical columns
    stacked_row_size = [792, 472]
    stacked_row_position = [10, 10]

    stacked_row_color_legend = {"display": True, "pos": [375, 1]}
    
    stacked_row_xaxis_properties = {"label": "Frequency", "textSize": 8, "tickfont": 8, "size": 55}
    stacked_row_yaxis_properties = {"label": "label1", "textSize": 8, "tickfont": 8, "size": 110}

    # Create stacked row plot
    stacked_row_plot = create_stacked_row_plot(
        stacked_row_title, stacked_row_params, stacked_row_size, stacked_row_position,
        stacked_row_color_legend, stacked_row_xaxis_properties, stacked_row_yaxis_properties
    )
    
    # Convert stacked row plot to JSON
    stacked_row_plot_json = convert_plot_to_json(stacked_row_plot)
    
    # SelectionDialogPlot parameters
    selection_dialog_plot_title = "Selection Dialog Plot Example"
    selection_dialog_params = ["param2", "param3", "param4", "param5", "param6"]  # Use any column
    selection_dialog_size = [792, 472]
    selection_dialog_position = [10, 500]

    # Create selection dialog plot
    selection_dialog_plot = create_selection_dialog_plot(
        selection_dialog_plot_title, selection_dialog_params, selection_dialog_size, selection_dialog_position
    )
    
    # Convert selection dialog plot to JSON
    selection_dialog_plot_json = convert_plot_to_json(selection_dialog_plot)

    view_config = {'initialCharts': {datasource_name: [stacked_row_plot_json, selection_dialog_plot_json]}}

    # Set views and serve project
    project.set_view(view_name, view_config)
    project.set_editable(True)
    # project.serve()

if __name__ == "__main__":
    main()