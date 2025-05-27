## This file creates a new MDV project named 'project' containing a view named 'default'.
## The datasource used is an h5ad file that was provided at runtime.
## The view 'default' shows one abundance box plot.
## The h5ad file is an AnnData object, and the `obs` attribute was used to extract the metadata for cells.
## The abundance box plot uses Disease on the x-axis and Abundance on the y-axis, with categorical variables Disease, Final_analysis, and Inflammation as parameters.
## This plot visualizes the distribution of sample IDs across disease states for different cell types, helping identify patterns specific to certain conditions.
## By examining the abundance of samples for each disease and cell type, one can detect shifts in cellular composition that correlate with disease pathology or severity. 
## For example, certain cell types may show increased abundance in diseased states compared to healthy ones, providing insights into disease-associated cellular changes.
## The plot also facilitates exploration of heterogeneity between patients, highlighting differences in cell type distributions that could suggest personalized treatment strategies or biological variability. 
## This visualization can be particularly useful for identifying cell types enriched in diseases like Crohn's disease, ulcerative colitis, or other inflammatory or autoimmune conditions.


import os
import json
import pandas as pd
import scanpy as sc
from mdvtools.mdvproject import MDVProject
from mdvtools.charts.abundance_box_plot import AbundanceBoxPlot
from mdvtools.charts.text_box_plot import TextBox

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

def create_text_box_plot(title, param, text, size, position, plot_id):
    """Create and configure a TextBoxPlot instance with the given parameters."""
    plot = TextBox(
        title=title,
        param=param,
        size=size,
        position=position,
        id=plot_id
    )
    
    plot.set_text(text) # setting the text
    
    return plot

def convert_plot_to_json(plot):
    """Convert plot data to JSON format."""
    return json.loads(json.dumps(plot.plot_data, indent=2).replace("\\", ""))

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
    
    # AbundanceBoxPlot parameters
    title = "Abundance box plot"
    params = ["Disease", "Sample ID", "Final_analysis"] #param1, param2 and param3 should all be categorical variables
    size = [615, 557]
    position = [341, 49]
    plot_id = "tGa0CF"
    x_axis_labels = [""]
    x_axis_title = ""
    y_axis_labels = ["Abundance"]
    y_axis_title = "Abundance"

    # TextBoxPlot parameters
    title_textbox = "Useful Information"
    param_textbox = []
    size_textbox = [615, 557]
    position_textbox = [341, 49]
    plot_id_textbox = "tGa0CF"
    text = "This plot visualizes the distribution of sample IDs across disease states for different cell types, helping identify patterns specific to certain conditions."
    
    # Create and configure plots
    abundance_box_plot = create_abundance_box_plot(
        title, params, size, position, plot_id, x_axis_labels, x_axis_title, y_axis_labels, y_axis_title
    )

    text_box_plot = create_text_box_plot(
        title_textbox, param_textbox, text, size_textbox, position_textbox, plot_id_textbox
    )
    
    # Convert plot to JSON and set view
    abundance_chart_json = convert_plot_to_json(abundance_box_plot)
    text_chart_json = convert_plot_to_json(text_box_plot)
    plots_view = {'initialCharts': {datasource_name: [text_chart_json]}}
    
    project.set_view(view_name, plots_view)
    project.set_editable(True)
    project.serve()

if __name__ == "__main__":
    main()


