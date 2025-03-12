import os
import json
import pandas as pd
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

def load_data(path):
    """Load data from the specified CSV file."""
    return pd.read_csv(path, low_memory=False)

def convert_plot_to_json(plot):
    """Convert plot data to JSON format."""
    return json.loads(json.dumps(plot.plot_data, indent=2).replace("\\\\", ""))

def main():
    """Main function to create the project and serve it."""
    # Constants
    project_path = os.path.expanduser('~/mdv/project')
    data_path = "path_to_data"
    view_name = "default"
    
    # Create project
    project = MDVProject(project_path, delete_existing=True)
    
    # Load data
    data_frame = load_data(data_path)

    # Set the correct data type to the "leiden" data source (imports as integer but it should be str to appear as a category)
    data_frame['leiden'] = data_frame['leiden'].apply(str)
    
    # Add datasource
    project.add_datasource(data_path, data_frame)
    
    # StackedRowChart parameters
    title = "Stacked Row Plot Example"
    params = ["leiden", "leiden"]
    size = [792, 472]
    position = [10, 10]

    #bandwith = 0.1
    #intervals = 40
    #scale = "0.001"

    legend_display = True
    #legend_position = [375,1]
              
    
    xaxis_properties = {"label": "n_genes_by_counts", 
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
    stackedrowchart_view = {'initialCharts': {data_path: [stacked_row_plot_json]}}
    
    project.set_view(view_name, stackedrowchart_view)
    project.set_editable(True)
    project.serve()

if __name__ == "__main__":
    main()
