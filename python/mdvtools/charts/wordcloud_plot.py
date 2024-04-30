# wordcloud_plot.py

from mdvtools.charts.base_plot import BasePlot

class WordcloudPlot(BasePlot):
    def __init__(self, title, param, wordSize, size, position, id=None):
        super().__init__(title, "row_chart", [param], size, position, id)
        self.plot_data["wordSize"] = wordSize
        self.plot_data["wordcloud"] = True

    def set_axis_properties(self, axis, properties):
        self.plot_data["axis"] = properties

    # Any additional methods specific to wordclouds can be added here
