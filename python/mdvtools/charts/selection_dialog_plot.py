from mdvtools.charts.base_plot import BasePlot

class SelectionDialogPlot(BasePlot):
    def __init__(self, title, params, size, position, id=None):
        super().__init__(title, "selection_dialog", params, size, position, id)


    # Additional methods for color scaling, log scale, etc., can be added here