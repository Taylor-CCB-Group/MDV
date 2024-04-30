from mdvtools.charts.base_plot import BasePlot

class DotPlot(BasePlot):
    def __init__(self, title, params, size, position, id=None):
        super().__init__(title, "dot_plot", params, size, position, id)

    def set_axis_properties(self, axis, properties):
        if "axis" not in self.plot_data:
            self.plot_data["axis"] = {}
        self.plot_data["axis"][axis] = properties

    def set_color_scale(self, log_scale=False):
        self.plot_data["color_scale"] = {"log": log_scale}

    def set_color_legend(self, display, position):
        self.plot_data["color_legend"] = {
            "display": display,
            "pos": position
        }

    def set_fraction_legend(self, display, position):
        self.plot_data["fraction_legend"] = {
            "display": display,
            "pos": position
        }

    def set_position(self, position):
        self.position = position

    # Additional methods for customization (e.g., tooltip visibility) can be added here
