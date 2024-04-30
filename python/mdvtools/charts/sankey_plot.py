from mdvtools.charts.base_plot import BasePlot

class SankeyPlot(BasePlot):
    def __init__(self, title, params, size, position, id=None):
        super().__init__(title, "sankey_chart", params, size, position, id)
        
    def set_axis_properties(self, axis, properties):
        self.plot_data[axis] = properties
