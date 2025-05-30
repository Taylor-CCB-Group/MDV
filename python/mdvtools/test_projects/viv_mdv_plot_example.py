# viv_mdv_plot_example.py

import os
import json
import pandas as pd
from mdvtools.mdvproject import MDVProject
from mdvtools.charts.viv_mdv_plot import VivMdvPlot

def create_viv_mdv_plot(title, params, size, position):
    """Create and configure a VivMdvPlot instance with the given parameters."""
    plot = VivMdvPlot(
        title=title,
        params=params,
        size=size,
        position=position
    )
    
    # Set additional configuration options
    plot.set_opacity(1)
    plot.set_tooltip(show=False)
    plot.set_region(title)
    plot.set_roi(0, 0, 5661, 9323)
    
    return plot

def load_data(path):
    """Load data from the specified CSV file."""
    return pd.read_csv(path, low_memory=False)

def convert_plot_to_json(plot):
    """Convert plot data to JSON format."""
    return json.loads(json.dumps(plot.plot_data, indent=2))

def main():
    """Main function to create the project and configure the plot."""
    # Define constants
    project_path = os.path.expanduser('~/mdv/project')
    data_path = '/Users/mariak/Documents/MDVmk/MDV/python/mdvtools/data/data_cells.csv'
    view_name = "default"
    
    # Create project
    project = MDVProject(project_path, view_name)
    
    # Create and configure the plot
    plot = create_viv_mdv_plot(
        title="Unvaccinated-Lesioned_SAMPLE_45_ROI_1",
        params=["x", "y", "leiden"],
        size=[1510, 767],
        position=[0, 0]
    )

    # Add plot to project
    project.add_datasource("cells", convert_plot_to_json(plot))

    # Convert project configuration to JSON and save it
    with open("viv_mdv_output.json", "w") as f:
        json.dump({
            "initialCharts": {"cells": [plot.to_dict()]},
            "dataSources": {"cells": {"layout": "gridstack", "panelWidth": 100}}
        }, f, indent=2)

if __name__ == "__main__":
    main()
