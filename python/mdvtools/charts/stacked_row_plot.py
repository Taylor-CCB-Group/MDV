from mdvtools.charts.base_plot import BasePlot

class StackedRowChart(BasePlot):
    def __init__(self, title, params, size, position, id=None):
        super().__init__(title, "stacked_row_chart", params, size, position, id)

    def set_axis_properties(self, axis, properties):
        if "axis" not in self.plot_data:
            self.plot_data["axis"] = {}
        self.plot_data["axis"][axis] = properties

    def set_color_legend(self, display):
        self.plot_data["color_legend"] = {
            "display": display
        }

    def set_position(self, position):
        self.position = position

    # Additional methods for customization can be added here
