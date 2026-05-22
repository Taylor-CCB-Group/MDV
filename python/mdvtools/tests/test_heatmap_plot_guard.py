import pytest

from mdvtools.charts.heatmap_plot import HeatmapPlot


def test_heatmap_plot_requires_at_least_two_params():
    with pytest.raises(ValueError, match="requires params"):
        HeatmapPlot(
            title="t",
            params=["only_category"],
            size=[100, 100],
            position=[0, 0],
        )

