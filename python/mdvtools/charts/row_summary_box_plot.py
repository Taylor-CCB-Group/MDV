from mdvtools.charts.base_plot import BasePlot

class RowSummaryBox(BasePlot):
    def __init__(self, title, param, size, position, id=None):
        super().__init__(title, "row_summary_box", param, size, position, id)