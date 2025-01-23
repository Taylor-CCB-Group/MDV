import os
import json
import pandas as pd
import scanpy as sc
from mdvtools.mdvproject import MDVProject
from mdvtools.charts.stacked_row_plot import StackedRowChart

def create_stacked_row_plot(title, params, size, position, legend_display, xaxis_properties, yaxis_properties):
    """Create and configure a StackedRowChart instance with the given parameters."""
    plot = StackedRowChart(
        title=title,
        params=params,
        size=size,
        position=position
    )

    plot.set_color_legend(legend_display)
    plot.set_axis_properties("x", xaxis_properties)
    plot.set_axis_properties("y", yaxis_properties)

    return plot

def convert_plot_to_json(plot):
    """Convert plot data to JSON format."""
    return json.loads(json.dumps(plot.plot_data, indent=2).replace("\\\\", ""))

def main():
    """Main function to create the project and serve it."""
    # Constants
    project_path = os.path.expanduser('~/mdv/project')
    data_path = "path_to_data"
    view_name = "default"
    datasource_name = "datasource_name"
    
    # Create project
    project = MDVProject(project_path, delete_existing=True)
    
    # Load data
    adata = sc.read_h5ad(data_path)
    data_frame = pd.DataFrame(adata.obs)

    # Add datasource
    project.add_datasource(datasource_name, data_frame)
    
    # StackedRowChart parameters
    title = "Stacked Row Plot Example"
    params = ["param1", "param2"] #param1 and param2 should all be categorical variables
    size = [792, 472]
    position = [10, 10]

    legend_display = True              
    
    xaxis_properties = {"label": "label1", 
             "textSize": 13, 
             "tickfont": 10
    }

    yaxis_properties = {"label": "density", 
             "textSize": 13, 
             "tickfont": 10
    }

    # Create plot
    stacked_row_plot = create_stacked_row_plot(title, params, size, position, legend_display, xaxis_properties, yaxis_properties)
    
    # Convert plot to JSON and set view
    stacked_row_plot_json = convert_plot_to_json(stacked_row_plot)
    stackedrowchart_view = {'initialCharts': {datasource_name: [stacked_row_plot_json]}}
    
    project.set_view(view_name, stackedrowchart_view)
    project.set_editable(True)
    project.serve()

if __name__ == "__main__":
    main()
