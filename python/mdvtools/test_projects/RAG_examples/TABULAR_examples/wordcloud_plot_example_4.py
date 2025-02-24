import os
import pandas as pd
import json
from mdvtools.mdvproject import MDVProject
from mdvtools.charts.wordcloud_plot import WordcloudPlot


def create_wordcloud_plot(title, param, wordSize, size, position, axis_settings):
    """Create and configure a WordcloudPlot instance with the given parameters."""
    plot = WordcloudPlot(
        title=title,
        param=param, #the param has to be just one categorical variable in the form of a string
        wordSize=wordSize,
        size=size,
        position=position
    )
    
    plot.set_axis_properties("x", axis_settings)  # axis settings
    
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
    
    # Wordcloud parameters
    title = "Wordcloud Plot Example"
    param = "param1" #param1 should be a categorical variable
    wordSize = 55  
    size = [750, 500]  
    position = [30, 30]
    
    axis_settings = {
        'textSize': 13,
        'label': "Axis label",
        'tickfont': 10
    }
    
    # Create plot
    wordcloud_plot = create_wordcloud_plot(title, param, wordSize, size, position, axis_settings)

    # Convert plot to JSON and set view
    wordcloud_plot_json = convert_plot_to_json(wordcloud_plot)
    wordcloud_plot_view = {'initialCharts': {datasource_name: [wordcloud_plot_json]}}
    
    project.set_view(view_name, wordcloud_plot_view)
    project.set_editable(True)
    project.serve()

if __name__ == "__main__":
    main()
