# density_scatter_plot.py

from mdvtools.charts.base_plot import BasePlot

class DensityScatterPlot(BasePlot):
    def __init__(self, title, params, size, position, id=None):
        super().__init__(title, "density_scatter_plot", params, size, position, id)

    def set_x_axis(self, label, size, textsize, tickfont):
        self.plot_data["x_axis"] = {
            "label": label,
            "size": size,
            "textSize": textsize,
            "tickFont": tickfont
        }

    def set_y_axis(self, label, size, textsize, tickfont, rotate_labels):
        self.plot_data["y_axis"] = {
            "label": label,
            "size": size,
            "textSize": textsize,
            "tickFont": tickfont,
            "rotateLabels": rotate_labels
        }


    def set_visuals(self, default_color, radius, opacity, log_color_scale, fallback_on_zero, background_color, bandwidth, intensity, opacity_cnt):
        self.plot_data.update({
            "defaultColor": default_color,
            "radius": radius,
            "opacity": opacity,
            "logColorScale": log_color_scale,
            "fallbackOnZero": fallback_on_zero,
            "backgroundColor": background_color,
            "contour_bandwidth": bandwidth,
            "contour_intensity": intensity,
            "contour_opacity": opacity_cnt
        })

    # Additional methods specific to density scatter plots can be added here
