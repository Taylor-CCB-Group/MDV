import os
import pandas as pd
import scanpy as sc
import sys
import  numpy as np
from mdvtools.mdvproject import MDVProject
from mdvtools.charts.heatmap_plot import HeatmapPlot
from mdvtools.charts.histogram_plot import HistogramPlot
from mdvtools.charts.dot_plot import DotPlot
from mdvtools.charts.box_plot import BoxPlot
from mdvtools.charts.scatter_plot_3D import ScatterPlot3D
from mdvtools.charts.row_chart import RowChart
from mdvtools.charts.scatter_plot import ScatterPlot
from mdvtools.charts.abundance_box_plot import AbundanceBoxPlot
from mdvtools.charts.stacked_row_plot import StackedRowChart
from mdvtools.charts.ring_chart import RingChart
from mdvtools.charts.violin_plot import ViolinPlot
from mdvtools.charts.multi_line_plot import MultiLinePlot
from mdvtools.charts.pie_chart import PieChart

import json 




def load_data(path):
    #Load data from the specified CSV file.
    return pd.read_csv(path, low_memory=False)

def convert_plot_to_json(plot):
    #Convert plot data to JSON format.
    return json.loads(json.dumps(plot.plot_data, indent=2).replace("\\", ""))
    

def create_bar_chart(title, params, size, position, color_legend, xaxis_properties, yaxis_properties):
    """Create and configure a StackedRowChart instance as a bar chart with the given parameters."""
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
    """Convert plot data to JSON format."""
    return json.loads(json.dumps(plot.plot_data, indent=2).replace("\\\\", ""))

def main():
    """Main function to create the project and serve it."""
    # Constants
    project_path = os.path.expanduser('~/mdv/project')
    view_name = "default"
    
    # Load data
    data_path = "file_path"
    adata = sc.read_h5ad(data_path)
    cells_df = pd.DataFrame(adata.obs)
    datasource_name = "datasource_name"
    
    # Create project
    project = MDVProject(project_path, delete_existing=True)
    
    # Add datasource
    project.add_datasource(datasource_name, cells_df)
    
    # Bar chart parameters
    title = "Number of samples per cell type"
    params = ["final_analysis", "sample_id"]  # Using 'sample_id' for samples and 'CellsLoaded' for the number of cells
    size = [792, 472]
    position = [10, 10]

    color_legend = {"display" : True,
                    "pos" : [375,1]}
    
    xaxis_properties = {"label": "Cell type", 
                        "textSize": 13, 
                        "tickfont": 10
    }

    yaxis_properties = {"label": "Percent of Cells", 
                        "textSize": 13, 
                        "tickfont": 10
    }

    # Create bar chart
    bar_chart = create_bar_chart(title, params, size, position, color_legend, xaxis_properties, yaxis_properties)
    
    # Convert plot to JSON and set view
    bar_chart_json = convert_plot_to_json(bar_chart)
    bar_chart_view = {'initialCharts': {datasource_name: [bar_chart_json]}}
    
    project.set_view(view_name, bar_chart_view)
    project.set_editable(True)
    project.serve()

if __name__ == "__main__":
    main()