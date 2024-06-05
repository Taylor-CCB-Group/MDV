# Assuming this should be placed in a new file, maybe scatter_plot.py

from mdvtools.charts.base_plot import BasePlot

class ScatterPlot3D(BasePlot):
    def __init__(self, title, params, size, position, id=None, default_color="#377eb8", brush="default", center=[0, 0, 0], on_filter="hide", radius=5, opacity=0.8, axis_scales=[1, 1, 1], camera=None):
        super().__init__(title, "wgl_3d_scatter_plot", params, size, position, id)
        self.plot_data["default_color"] = default_color
        self.plot_data["brush"] = brush
        self.plot_data["center"] = center
        self.plot_data["on_filter"] = on_filter
        self.plot_data["radius"] = radius
        self.plot_data["opacity"] = opacity
        self.plot_data["axis_scales"] = axis_scales
        self.plot_data["camera"] = camera if camera is not None else {"distance": 1000, "theta": 0, "phi": 0}

    # Set the color for the plot
    def set_default_color(self, color):
        self.plot_data["default_color"] = color

    # Set the brush type
    def set_brush(self, brush):
        self.plot_data["brush"] = brush

    # Set the center for the 3D plot
    def set_center(self, center):
        self.plot_data["center"] = center

    # Set the filter action
    def set_on_filter(self, on_filter):
        self.plot_data["on_filter"] = on_filter

    # Set the radius for the scatter points
    def set_radius(self, radius):
        self.plot_data["radius"] = radius

    # Set the opacity for the scatter points
    def set_opacity(self, opacity):
        self.plot_data["opacity"] = opacity

    # Set the scales for each axis
    def set_axis_scales(self, axis_scales):
        self.plot_data["axis_scales"] = axis_scales

    # Set the camera position for the 3D plot
    def set_camera(self, camera):
        self.plot_data["camera"] = camera

    def set_color_by(self, color_by):
        self.plot_data["color_by"] = color_by

    # Additional methods specific to 3D scatter plots can be added here
