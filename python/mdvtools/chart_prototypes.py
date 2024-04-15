"""
As of this writing, these are extracted from project-specific scripts
where they are used for generating views.
Expecting to develop into something with a very different API quite soon - this is highly provisional/unstable.
"""

import copy
import random
import string
from typing import Literal

selection_chart_proto = {
    "filters": {"sampleID": {"category": ["21455"], "operand": None}},
    "id": "590SPa",
    "legend": "",
    "param": ["sampleID"],
    "size": [305, 155],
    "title": "",
    "type": "selection_dialog",
    "position": [0, 766],
    "gsposition": [0, 5],
    "gssize": [6, 1],
}

stat_chart_proto = {
    "axis": {
        "x": {
            "label": "xMin",
            "size": 30,
            "textSize": 13,
            "textsize": 13,
            "tickfont": 10,
        },
        "y": {
            "label": "yMin",
            "size": 45,
            "textSize": 13,
            "textsize": 13,
            "tickfont": 10,
        },
    },
    "brush": "poly",
    "color_by": "gmax_Periostin-Neutrophil",
    "color_legend": {"display": True, "pos": [45, 10]},
    "default_color": "#377eb8",
    "id": "2JxLCo",
    "legend": "",
    "log_color_scale": True,
    "on_filter": "hide",
    "opacity": 0.8,
    "param": ["xMin", "yMin"],
    "radius": 5.12,
    "size": [203, 311],
    "title": "feature3",
    "tooltip": {"show": True},
    "trim_color_scale": "0.05",
    "type": "wgl_scatter_plot",
    "position": [0, 0],
    "gsposition": [8, 0],
    "gssize": [12, 5],
}

violin_proto = {
    "axis": {
        "x": {
            "label": "Analysis type",
            "size": 30,
            "textSize": 13,
            "textsize": 13,
            "tickfont": 10,
        },
        "y": {
            "label": "importance",
            "size": 45,
            "textSize": 13,
            "textsize": 13,
            "tickfont": 10,
        },
    },
    "band_width": 0.00038430866,
    "brush": "poly",
    "color_by": "feature",
    "color_legend": {"display": False, "pos": [45, 10]},
    "default_color": "#377eb8",
    "gsposition": [0, 3],
    "gssize": [12, 2],
    "id": "lu3qYD",
    "intervals": 20,
    "legend": "",
    "on_filter": "hide",
    "opacity": 0.8,
    "param": ["Analysis type", "importance"],
    "radius": 5,
    "size": [390, 291],
    "title": "Analysis type x importance",
    "tooltip": {"column": "feature", "show": True},
    "type": "violin_plot",
}

feature_import_proto = {
    "title": "calculated importance of feature for patient",
    "legend": "",
    "type": "table_chart",
    "param": [
        "feature",
        "importance",
    ],
    "include_index": False,
    "id": "le549N",
    "size": [1074, 980],
    "column_widths": {"feature": 386},
    "gsposition": [0, 0],
    "gssize": [12, 3],
}

contribution_heatmap_proto = {
    "title": "13309-cell_pair_contribution",
    "legend": "",
    "type": "single_heat_map",
    "param": ["Cell Type 1", "Cell Type 2", "cell_pair_contribution", "sampleID"],
    "pivot": "13309",
    "value": "cell_pair_contribution",
    "category": "13309",
    "id": "b1eIsa",
    "size": [300, 300],
    "axis": {
        "x": {"textSize": 13, "label": "", "size": 25, "tickfont": 10},
        "y": {"textSize": 13, "label": "", "size": 40, "tickfont": 10},
    },
    "color_scale": {"log": False, "min_max": [0, 0.3117565488815312]},
    "color_legend": {"display": False, "pos": [48, 164]},
    "position": [10, 10],
    "gsposition": [0, 0],
    "gssize": [12, 3],
}

viv_proto = {
    "radius": 6,
    "opacity": 1,
    "color_by": "Celltype",
    "color_legend": {"display": True, "pos": [3, 48]},
    "tooltip": {"show": True, "column": "Celltype"},
    "background_filter": {"category": "48639_immune", "column": "panelID"},
    "channel": 0,
    "gsposition": [0, 0],
    "gssize": [12, 3],
    "id": "mLfaWl",
    "legend": "",
    "param": ["x", "y", "Celltype"],
    "region": "48639_immune",
    "roi": {"max_x": 1262, "max_y": 1548, "min_x": 0, "min_y": 0},
    "size": [445, 483],
    "title": "48639_immune",
    "type": "VivMdvRegionReact",
    "viv": {
        "channels": [{"name": "DAPI"}],
        "image_properties": {
            "channelsVisible": [],
            "colors": [],
            "contrastLimits": [],
            "domains": [],
            "ids": [],
            "image": 0,
            "loader": [{"labels": [], "shape": []}],
            "selections": [],
        },
    },
}

stat_viv_proto = {
    "radius": 6,
    "opacity": 0.25,
    "color_by": "gmax_Periostin-Neutrophil",
    "color_legend": {"display": True, "pos": [3, 48]},
    "tooltip": {"show": True, "column": "gmax_Periostin-Neutrophil"},
    "background_filter": {"category": "48639", "column": "sampleID"},
    "channel": 0,
    "gsposition": [0, 0],
    "gssize": [12, 5],
    "id": "mLfaWl",
    "legend": "",
    "param": ["x", "y", "gmax_Periostin-Neutrophil"],
    "region": "48639",
    "roi": {"max_x": 1262, "max_y": 1548, "min_x": 0, "min_y": 0},
    "size": [445, 483],
    "title": "48639",
    "type": "VivMdvRegionReact",
    "viv": {
        "channels": [{"name": "DAPI"}],
        "image_properties": {
            "channelsVisible": [],
            "colors": [],
            "contrastLimits": [],
            "domains": [],
            "ids": [],
            "image": 0,
            "loader": [{"labels": [], "shape": []}],
            "selections": [],
        },
    },
    "trim_color_scale": "0.05",
}

multi_line_proto = {
    "title": "CD103",
    "legend": "",
    "type": "multi_line_chart",
    "param": ["CD103", "sample"],
    "axis": {
        "x": {
            "size": 30,
            "label": "CD103",
            "textsize": 13,
            "textSize": 13,
            "tickfont": 10,
        },
        "y": {
            "size": 45,
            "label": "density",
            "textsize": 13,
            "textSize": 13,
            "tickfont": 10,
        },
    },
    "id": "NjoSX8",
    "size": [250, 309],
    "color_legend": {"display": False, "pos": [456, 14]},
    "band_width": 0.1,
    "intervals": 100,
    "stacked": False,
    "fill": False,
    "scaletrim": "0.05",
    "gsposition": [1, 5],
    "gssize": [1, 2],
}


chart_prototypes = {
    "selection": selection_chart_proto,
    "stats": stat_chart_proto,
    "stats_viv": stat_viv_proto,
    "feature": feature_import_proto,
    "contribution_heatmap": contribution_heatmap_proto,
    "viv": viv_proto,
    "violin": violin_proto,
    "multi_line": multi_line_proto,
}

# for hinting, may not maintain this, probably use a different general approach later (more class-based etc)
ChartKeys = Literal[
    "selection",
    "stats",
    "stats_viv",
    "feature",
    "contribution_heatmap",
    "viv",
    "violin",
    "multi_line",
]


def make_chart(proto_name: ChartKeys):
    c = copy.deepcopy(chart_prototypes[proto_name])
    c["id"] = str("".join(random.choices(string.ascii_letters, k=6)))
    return c


def selection_chart(sampleID, column="sampleID"):
    selection_chart = make_chart("selection")
    # selection_chart['filters']['sampleID']['category'][0] = sampleID
    selection_chart["filters"] = {column: {"category": [sampleID], "operand": None}}
    selection_chart["param"] = [column]
    return selection_chart


def image_chart(imageID, region, color_by, image_column="panelID"):
    c = make_chart("viv")
    c["title"] = c["region"] = c["background_filter"]["category"] = imageID
    c["roi"]["max_x"] = region["width"]
    c["roi"]["max_y"] = region["height"]
    c["param"][2] = c["color_by"] = c["tooltip"]["column"] = color_by
    c["background_filter"]["column"] = image_column
    return c
