import os
import json
import pandas as pd
import scanpy as sc
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

def convert_plot_to_json(plot):
    """Convert plot data to JSON format."""
    return json.loads(json.dumps(plot.plot_data, indent=2).replace("\\", ""))

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
    
    # AbundanceBoxPlot parameters
    title = "Example title"
    params = ["param1", "param2", "param3"] #param1, param2 and param3 should all be categorical variables. Params take only 3 arguments.
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
    abundance_view = {'initialCharts': {datasource_name: [abundance_chart_json]}}
    
    project.set_view(view_name, abundance_view)
    project.set_editable(True)
    project.serve()

if __name__ == "__main__":
    main()


