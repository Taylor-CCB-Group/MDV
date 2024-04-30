from mdvtools.charts.base_plot import BasePlot

class HistogramPlot(BasePlot):
    def __init__(self, title, param, bin_number, display_min, display_max, size, position, id=None):
        super().__init__(title, "bar_chart", param, size, position, id)
        self.plot_data["bin_number"] = bin_number
        self.plot_data["display_min"] = display_min
        self.plot_data["display_max"] = display_max

    def set_x_axis(self, size, label, textsize, tickfont):
        self.plot_data["x_axis"] = {
            "size": size,
            "label": label,
            "textSize": textsize,
            "tickfont": tickfont
        }

    def set_y_axis(self, size, label, textsize, tickfont, rotate_labels):
        self.plot_data["y_axis"] = {
            "size": size,
            "label": label,
            "textSize": textsize,
            "tickfont": tickfont,
            "rotate_labels": rotate_labels
        }
