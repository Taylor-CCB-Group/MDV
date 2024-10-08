import os
import json
import pandas as pd
from mdvtools.mdvproject import MDVProject
from mdvtools.charts.abundance_box_plot import AbundanceBoxPlot

def create_abundance_box_plot(title, params, size, position, plot_id, x_axis_labels, x_axis_title, y_axis_labels, y_axis_title):
    """Create and configure an AbundanceBoxPlot instance with the given parameters."""
    plot = AbundanceBoxPlot(
        title=title,
        params=params,
        size=size,
        position=position,
        id=plot_id
    )
    
    plot.set_x_axis(x_axis_labels, x_axis_title)  # x-axis labels and title
    plot.set_y_axis(y_axis_labels, y_axis_title)  # y-axis labels and title
    
    return plot

def load_data(path):
    """Load data from the specified CSV file."""
    return pd.read_csv(path, low_memory=False)

def convert_plot_to_json(plot):
    """Convert plot data to JSON format."""
    return json.loads(json.dumps(plot.plot_data, indent=2).replace("\\", ""))

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
    
    # Add datasource
    project.add_datasource(data_path, data_frame)
    
    # AbundanceBoxPlot parameters
    title = "PBMC3k gene expression"
    params = ["mt", "mt", "mt"]
    size = [615, 557]
    position = [341, 49]
    plot_id = "tGa0CF"
    x_axis_labels = [""]
    x_axis_title = ""
    y_axis_labels = ["Abundance"]
    y_axis_title = "Abundance"
    
    # Create and configure plot
    abundance_box_plot = create_abundance_box_plot(
        title, params, size, position, plot_id, x_axis_labels, x_axis_title, y_axis_labels, y_axis_title
    )
    
    # Convert plot to JSON and set view
    abundance_chart_json = convert_plot_to_json(abundance_box_plot)
    abundance_view = {'initialCharts': {data_path: [abundance_chart_json]}}
    
    project.set_view(view_name, abundance_view)
    project.set_editable(True)
    project.serve()

if __name__ == "__main__":
    main()


