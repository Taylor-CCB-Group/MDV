# scatter_plot.py

from mdvtools.charts.base_plot import BasePlot

class ScatterPlot(BasePlot):
    def __init__(self, title, params, size, position, id=None):
        super().__init__(title, "wgl_scatter_plot", params, size, position, id)

    def set_default_color(self, color):
        self.plot_data["default_color"] = color

    def set_brush(self, brush_type):
        self.plot_data["brush"] = brush_type

    def set_filter(self, filter_type):
        self.plot_data["on_filter"] = filter_type

    def set_radius(self, radius):
        self.plot_data["radius"] = radius

    def set_opacity(self, opacity):
        self.plot_data["opacity"] = opacity

    def set_color_by(self, color_by):
        self.plot_data["color_by"] = color_by

    def set_color_legend(self, display, position):
        self.plot_data["color_legend"] = {
            "display": display,
            "pos": position
        }

    def set_tooltip(self, show):
        self.plot_data["tooltip"] = {"show": show}
    
    # Any additional methods specific to scatter plots can be added here

