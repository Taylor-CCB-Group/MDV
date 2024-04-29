from mdvtools.charts.base_plot import BasePlot

class TablePlot(BasePlot):
    def __init__(self, title, params, size, position, id=None, index=False,):
        super().__init__(title, "table_chart", params, size, position, id)
        self.plot_data["include_index"] = index


    # Any additional methods specific to wordclouds can be added here
