import os
import json
import pandas as pd
import scanpy as sc
from mdvtools.mdvproject import MDVProject
from mdvtools.charts.ring_chart import RingChart

def create_ring_chart(title, param, size, position):
    """Create and configure a RingChart instance with the given parameters."""
    plot = RingChart(
        title=title,
        param=param,
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
    view_name = "default"
    datasource_name = "datasource_name"
    
    # Create project
    project = MDVProject(project_path, delete_existing=True)
    
    # Load data
    adata = sc.read_h5ad(data_path)
    data_frame = pd.DataFrame(adata.obs)
    
    # Add datasource
    project.add_datasource(datasource_name, data_frame)
    
    # Ring chart parameters
    title = "Ring Chart Example"
    param = "param1" #param1 should be categorical
    size = [792, 472]
    position = [10, 10]
    
    # Create plot
    plot = create_ring_chart(title, param, size, position)
    
    # Convert plot to JSON and set view
    ring_chart_json = convert_plot_to_json(plot)
    Ringchart_view = {'initialCharts': {datasource_name: [ring_chart_json]}}
    
    project.set_view(view_name, Ringchart_view)
    project.set_editable(True)
    project.serve()

if __name__ == "__main__":
    main()
