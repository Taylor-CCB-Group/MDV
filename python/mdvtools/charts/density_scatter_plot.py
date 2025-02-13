from mdvtools.charts.base_plot import BasePlot


class DensityScatterPlot(BasePlot):
    def __init__(self, title, param, size, position, id=None):
        super().__init__(title, "density_scatter_plot", param, size, position, id)

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
        self.plot_data["color_legend"] = {"display": display, "pos": position}

    def set_tooltip(self, show):
        self.plot_data["tooltip"] = {"show": show}

    def set_category1(self, category1):
        self.plot_data["category1"] = category1

    def set_category2(self, category2):
        self.plot_data["category2"] = category2

    # Any additional methods specific to scatter plots can be added here
