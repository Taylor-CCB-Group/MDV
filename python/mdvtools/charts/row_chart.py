from mdvtools.charts.base_plot import BasePlot


class RowChart(BasePlot):
    def __init__(self, title, param, size, position, id=None, **kwargs):
        super().__init__(title, "row_chart", [param], size, position, id)

    def set_axis_properties(self, axis, properties):
        if "axis" not in self.plot_data:
            self.plot_data["axis"] = {}
        self.plot_data["axis"][axis] = properties

    # Additional methods for customization can be added here
