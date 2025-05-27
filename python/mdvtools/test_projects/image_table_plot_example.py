import os
import json
from mdvtools.mdvproject import MDVProject
from mdvtools.charts.image_table_plot import ImageTableChart


def create_image_table_chart(title, param, size, position, image_settings, grid_settings, legend=""):
    """Create and configure an ImageTableChart instance with the given parameters."""
    chart = ImageTableChart(
        title=title,
        param=param,
        size=size,
        position=position,
        #region= region,
        id="Uzvm8X"
    )
    
    chart.set_region(region = "image_name")

    chart.set_images(
        base_url=image_settings["base_url"],
        image_key=image_settings["key_column"],
        image_type=image_settings["type"],
        image_name=image_settings["name"]
    )
    chart.set_image_width(image_settings.get("width", 100))
    chart.set_margins(**image_settings.get("margins", {"top_bottom": 10, "left_right": 10}))
    chart.set_image_label(image_settings.get("label", "DAPI"))
    
    chart.set_grid_position(grid_settings["position"])
    chart.set_grid_size(grid_settings["size"])
    chart.set_legend(legend)
    
    return chart

def convert_chart_to_json(chart):
    """Convert chart data to JSON format."""
    return json.loads(json.dumps(chart.plot_data, indent=2))

def main():
    """Main function to create the project and serve it."""
    # Constants
    project_path = os.path.expanduser('~/mdv/project')
    view_name = "default"
    
    # Create project
    project = MDVProject(project_path, delete_existing=True)

    # Chart parameters
    title = "test_rgb"
    param = ["x", "y", "leiden"]
    size = [1510, 767]
    position = [0, 0]
    region = "image_name"

    image_settings = {
        "base_url": "images/EXAMPLE",
        "key_column": "sample_id",
        "type": "VivMdvRegionReact",
        "name": "Users/mariak/Downloads/test_rgb.ome.tiff",
        "width": 100,
        "margins": {"top_bottom": 10, "left_right": 10},
        "label": "DAPI"
    }
    
    grid_settings = {
        "position": [0, 0],
        "size": [12, 5]
    }

    # Create chart
    image_table_chart = create_image_table_chart(title, param, size, position, region, image_settings, grid_settings)

    # Convert chart to JSON and set view
    image_table_chart_json = convert_chart_to_json(image_table_chart)
    chart_view = {'initialCharts': {"cells": [image_table_chart_json]}}
    
    project.set_view(view_name, chart_view)
    project.set_editable(True)
    project.serve()

if __name__ == "__main__":
    main()
