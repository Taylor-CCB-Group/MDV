from mdvtools.charts.base_plot import BasePlot

class MultiLinePlot(BasePlot):
    def __init__(self, title, params, size, position, id=None):
        super().__init__(title, "multi_line_chart", params, size, position, id)
        self.plot_data["band_width"] = 0.1  # Default value, can be set by a method
        self.plot_data["intervals"] = 40  # Default value, can be set by a method
        self.plot_data["scaletrim"] = "0.001"  # Default value, can be set by a method

    def set_x_axis(self, label, size, textsize, tickfont):
        self.plot_data["x_axis"] = {
            "label": label,
            "size": size,
            "textSize": textsize,
            "tickFont": tickfont
        }

    def set_y_axis(self, label, size, textsize, tickfont):
        self.plot_data["y_axis"] = {
            "label": label,
            "size": size,
            "textSize": textsize,
            "tickFont": tickfont
        }

    def set_bandwidth(self, band_width):
        self.plot_data["band_width"] = band_width

    def set_intervals(self, intervals):
        self.plot_data["intervals"] = intervals

    def set_scaletrim(self, scaletrim):
        self.plot_data["scaletrim"] = scaletrim

    def set_color_legend(self, display, pos):
        self.plot_data["color_legend"] = {
            "display": display,
            "pos": pos
        }
