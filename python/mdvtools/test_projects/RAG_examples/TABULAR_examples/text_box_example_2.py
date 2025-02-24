import os
import json
import pandas as pd
from mdvtools.mdvproject import MDVProject
from mdvtools.charts.text_box_plot import TextBox

def create_text_box_plot(title, param, text, size, position):
    """Create and configure a TextBox instance with the given parameters."""
    plot = TextBox(
        title=title,
        param=param,
        size=size,
        position=position
    )

    plot.set_text(text)  # set the text

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
    datasource_name = "datasource_name"
    
    # Create project
    project = MDVProject(project_path, delete_existing=True)
    
    # Load data
    data_frame = load_data(data_path)
    
    # Add datasource
    project.add_datasource(datasource_name, data_frame)
    
    # TextBox parameters
    title = "Example title"
    param = [] #param should be empty.
    text = "This is a text box. Protein interactions reveal key molecular mechanisms."
    size = [850, 550]
    position = [40, 40]
   
    # Create plot
    plot = create_text_box_plot(title, param, text, size, position)
    
    # Convert plot to JSON and set view
    textboxplot_chart_json = convert_plot_to_json(plot)
    textboxplot_view = {'initialCharts': {datasource_name: [textboxplot_chart_json]}}
    
    project.set_view(view_name, textboxplot_view)
    project.set_editable(True)
    project.serve()

if __name__ == "__main__":
    main()
