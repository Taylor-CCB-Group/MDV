from mdvtools.charts.base_plot import BasePlot

class ImageTableChart(BasePlot):

    def __init__(self, title, params, size, position, id=None):
        super().__init__(title, "image_table_chart", params, size, position, id)

    def set_images(self, base_url, image_key, image_type):
        self.plot_data["images"] = {
            "base_url": base_url, #base_url is the path to the images the method 'add_image_set' is adding to the project, and it is always images/EXAMPLE
            "key_column": image_key,
            "type": image_type
        }

    def set_image_width(self, width):
        self.plot_data["image_width"] = width

    def set_margins(self, top_bottom, left_right):
        self.plot_data["margins"] = {
            "top_bottom": top_bottom,
            "left_right": left_right
        }

    def set_image_label(self, label):
        self.plot_data["image_label"] = label

    def set_grid_position(self, gsposition):
        self.plot_data["gsposition"] = gsposition

    def set_grid_size(self, gssize):
        self.plot_data["gssize"] = gssize

    def set_legend(self, legend):
        self.plot_data["legend"] = legend