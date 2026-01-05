# heatmap_plot.py

from mdvtools.charts.base_plot import BasePlot


class HeatmapPlot(BasePlot):
    def __init__(self, title, params, size, position, id=None, **kwargs):
        super().__init__(title, "heat_map", params, size, position, id, **kwargs)

    def set_method(self, method):
        self.plot_data["method"] = method

    def set_color_scale(self, color_scale):
        self.plot_data["color_scale"] = color_scale

    def set_intensity_scale(self, intensity_scale):
        self.plot_data["intensity_scale"] = intensity_scale

    def set_opacity(self, opacity):
        self.plot_data["opacity"] = opacity

    def set_color_range(self, color_range):
        self.plot_data["color_range"] = color_range

    def set_tooltip(self, show):
        self.plot_data["tooltip"] = {"show": show}

    def set_x_axis(self, axis_labels, axis_title):
        self.plot_data["x_axis"] = {"labels": axis_labels, "title": axis_title}

    def set_y_axis(self, axis_labels, axis_title):
        self.plot_data["y_axis"] = {"labels": axis_labels, "title": axis_title}

    def set_axis(self, xtextSize, ytextSize, xsize, ysize, xtickfont, ytickfont, xrotate_status, yrotate_status):
        self.plot_data["axis"] = {"x": {"textSize": xtextSize, "size": xsize, "tickfont": xtickfont, "rotate_labels": xrotate_status}, "y": {"textSize": ytextSize, "size": ysize, "tickfont": ytickfont, "rotate_labels": yrotate_status}}

    # Any additional methods specific to heatmaps can be added here
