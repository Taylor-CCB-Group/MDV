# violin_plot.py

from mdvtools.charts.base_plot import BasePlot

class ViolinPlot(BasePlot):
    def __init__(self, title, params, size, position, id=None):
        super().__init__(title, "violin_plot", params, size, position, id)

    def set_points(self, points):
        self.plot_data["points"] = points

    def set_box_visible(self, visible):
        self.plot_data["box_visible"] = visible

    def set_line_color(self, color):
        self.plot_data["line_color"] = color

    def set_fill_color(self, color):
        self.plot_data["fill_color"] = color

    def set_bandwidth(self, bandwidth):
        self.plot_data["bandwidth"] = bandwidth

    def set_meanline_visible(self, visible):
        self.plot_data["meanline_visible"] = visible

    def set_color_scale(self, scale):
        self.plot_data["color_scale"] = scale

    def set_tooltip(self, show, content=None):
        self.plot_data["tooltip"] = {
            "show": show,
            "content": content if content else "default"
        }

    # Additional methods specific to violin plots can be added here
