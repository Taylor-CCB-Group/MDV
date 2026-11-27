from mdvtools.template_views.create_question_views import spread_layout


def box(x: int, y: int, w: int, h: int) -> dict:
    return {"position": [x, y], "size": [w, h]}


def edges(chart: dict) -> tuple[int, int, int, int]:
    (x, y), (w, h) = chart["position"], chart["size"]
    return x, y, x + w, y + h


def test_gaps_widen_rows_stay_aligned_and_width_is_kept():
    charts = [
        box(0, 0, 1000, 100),
        box(0, 105, 495, 300),
        box(500, 105, 500, 300),
        box(0, 410, 300, 200),
        box(305, 410, 695, 200),
    ]
    spread_layout({"initialCharts": {"cells": charts}}, gap=25)
    text, left, right, low_left, low_right = (edges(chart) for chart in charts)

    assert text[2] == right[2] == low_right[2] == 1000
    assert left[0] == low_left[0] == 0
    assert right[0] - left[2] in (24, 25, 26)
    assert low_right[0] - low_left[2] in (24, 25, 26)
    assert left[1] - text[3] == 25
    assert low_left[1] - left[3] == 25
    assert [chart["size"][1] for chart in charts] == [100, 300, 300, 200, 200]


def test_single_column_panel_only_grows_taller():
    charts = [box(0, 0, 270, 400), box(0, 405, 270, 455)]
    spread_layout({"initialCharts": {"genes": charts}}, gap=25)

    assert charts == [box(0, 0, 270, 400), box(0, 425, 270, 455)]
