# volcano_plot.py

from mdvtools.charts.base_plot import BasePlot


class VolcanoPlot(BasePlot):
    def __init__(self, title, params, size, position, id=None, **kwargs):
        super().__init__(
            title, "volcano_plot", params, size, position, id, **kwargs
        )

    def set_fc_threshold(self, threshold):
        """Set the fold change threshold (default: 1.0)"""
        self.plot_data["fc_threshold"] = threshold

    def set_pval_threshold(self, threshold):
        """Set the p-value threshold (default: 0.05)"""
        self.plot_data["pval_threshold"] = threshold

    def set_show_threshold_lines(self, show):
        """Show or hide threshold lines (default: True)"""
        self.plot_data["show_threshold_lines"] = show

    def set_volcano_color_mode(self, mode):
        """Set color mode: 'significance' or 'single' (default: 'significance')"""
        self.plot_data["volcano_color_mode"] = mode

    def set_up_color(self, color):
        """Set color for up-regulated points (default: '#d62728')"""
        self.plot_data["up_color"] = color

    def set_down_color(self, color):
        """Set color for down-regulated points (default: '#1f77b4')"""
        self.plot_data["down_color"] = color

    def set_ns_color(self, color):
        """Set color for non-significant points (default: '#7f7f7f')"""
        self.plot_data["ns_color"] = color

    def set_default_color(self, color):
        """Set default color for single color mode"""
        self.plot_data["default_color"] = color

    def set_brush(self, brush_type):
        """Set brush type: 'poly' or 'default'"""
        self.plot_data["brush"] = brush_type

    def set_filter(self, filter_type):
        """Set filter behavior"""
        self.plot_data["on_filter"] = filter_type

    def set_radius(self, radius):
        """Set point radius"""
        self.plot_data["radius"] = radius

    def set_opacity(self, opacity):
        """Set point opacity"""
        self.plot_data["opacity"] = opacity

    def set_tooltip(self, show, column=None):
        """Show or hide tooltips and optionally set the column to display"""
        tooltip_config = {"show": show}
        if column:
            tooltip_config["column"] = column
        self.plot_data["tooltip"] = tooltip_config

    def set_tooltip_column(self, column):
        """Set which column to display in tooltips (also used for labels)"""
        if "tooltip" not in self.plot_data:
            self.plot_data["tooltip"] = {"show": False}
        self.plot_data["tooltip"]["column"] = column

    def set_show_labels(self, show):
        """Show or hide gene labels on the plot (default: True)"""
        self.plot_data["show_labels"] = show

    def set_label_column(self, column):
        """Set which column to use for gene labels (independent of tooltip)"""
        self.plot_data["label_column"] = column

    def set_label_count(self, count):
        """Set number of top genes to label (default: 10)"""
        self.plot_data["label_count"] = count
