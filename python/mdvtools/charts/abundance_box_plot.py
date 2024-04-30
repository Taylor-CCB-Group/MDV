from mdvtools.charts.base_plot import BasePlot

class AbundanceBoxPlot(BasePlot):
    def __init__(self, title, params, size, position, id=None):
        super().__init__(title, "custom_box_plot", params, size, position, id)

    def set_grouping(self, grouping):
        self.plot_data["grouping"] = grouping

    def set_log_transform(self, log_transform):
        self.plot_data["log_transform"] = log_transform

    def set_box_color(self, box_color):
        self.plot_data["box_color"] = box_color

    def set_whisker_color(self, whisker_color):
        self.plot_data["whisker_color"] = whisker_color

    def set_median_color(self, median_color):
        self.plot_data["median_color"] = median_color

    def set_outlier_color(self, outlier_color):
        self.plot_data["outlier_color"] = outlier_color

    def set_tooltip(self, show):
        self.plot_data["tooltip"] = {"show": show}

    def set_x_axis(self, axis_labels, axis_title):
        self.plot_data["x_axis"] = {
            "labels": axis_labels,
            "title": axis_title
        }

    def set_y_axis(self, axis_labels, axis_title):
        self.plot_data["y_axis"] = {
            "labels": axis_labels,
            "title": axis_title
        }

    # Any additional methods specific to abundance box plots can be added here
