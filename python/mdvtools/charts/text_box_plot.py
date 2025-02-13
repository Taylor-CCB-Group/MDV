from mdvtools.charts.base_plot import BasePlot

class TextBox(BasePlot):
    def __init__(self, title, param, size, position, id=None):
        super().__init__(title, "text_box_chart", param, size, position, id)

    def set_text(self, text):
        self.plot_data["text"] = text
