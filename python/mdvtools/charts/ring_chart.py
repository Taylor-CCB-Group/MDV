# ring_chart.py

from mdvtools.charts.base_plot import BasePlot

class RingChart(BasePlot):
    def __init__(self, title, param, size, position, id=None):
        super().__init__(title, "ring_chart", [param], size, position, id)
        self.plot_data["axis"] = {"x": {}}

    def set_legend(self, legend):
        self.plot_data["legend"] = legend

    def set_param(self, param):
        self.plot_data["param"] = param

    def set_axis_properties(self, axis, properties):
        if "axis" not in self.plot_data:
            self.plot_data["axis"] = {}
        self.plot_data["axis"][axis] = properties

    def set_text_size(self, text_size):
        self.plot_data["axis"]["x"]["textSize"] = text_size

    def set_label(self, label):
        self.plot_data["axis"]["x"]["label"] = label

    def set_axis_size(self, axis_size):
        self.plot_data["axis"]["x"]["size"] = axis_size

    def set_tick_font(self, tick_font):
        self.plot_data["axis"]["x"]["tickfont"] = tick_font

    #def set_hole_size(self, hole_size):
    #    self.hole_size = hole_size

    
    # Additional methods specific to ring charts can be added here
