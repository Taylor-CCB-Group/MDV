import os
import pandas as pd
import scanpy as sc
from mdvtools.mdvproject import MDVProject
from mdvtools.charts.box_plot import BoxPlot
import json
    

def create_box_plot(title, params, size, position, plot_id):
    """Create and configure a BoxPlot instance with the given parameters."""
    plot = BoxPlot(
        title=title,
        params=params,
        size=size,
        position=position,
        id=plot_id
    )
    return plot

def convert_plot_to_json(plot):
    """Convert plot data to JSON format."""
    return json.loads(json.dumps(plot.plot_data, indent=2).replace("\\", ""))

def main():
    """Main function to create the project and serve it."""
    # Constants
    project_path = os.path.expanduser('~/mdv/project')
    view_name = "Disease vs Total Counts Box Plot"
    
    # Load data
    data_path = "file_path"
    adata = sc.read_h5ad(data_path)
    cells_df = pd.DataFrame(adata.obs)
    cells_df.name = 'cells'
    
    # Create project
    project = MDVProject(project_path, delete_existing=True)
    
    # Add datasource
    project.add_datasource('cells', cells_df)
    
    # BoxPlot parameters
    box_title = "Total Counts per Disease"
    box_params = ["Disease", "total_counts"]
    box_size = [615, 557]
    box_position = [50, 50]
    box_plot_id = "boxPlot1"
    
    # Create box plot
    box_plot = create_box_plot(
        box_title, box_params, box_size, box_position, box_plot_id
    )
    
    # Convert plot to JSON and set view
    box_plot_json = convert_plot_to_json(box_plot)
    view_config = {'initialCharts': {'cells': [box_plot_json]}}
    
    project.set_view(view_name, view_config)
    project.set_editable(True)
    project.serve()

if __name__ == "__main__":
    main()