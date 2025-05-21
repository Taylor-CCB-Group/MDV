from mdvtools.mdvproject import MDVProject
import pandas as pd
import numpy as np
import random
import string
import os

rng = np.random.default_rng()
nrm = rng.normal

n = int(1e6)

x = nrm(size=n)
y = nrm(size=n)
z = nrm(size=n)

# generate some random text data
text = np.random.choice(["a", "b", "c"], n)
choices = ["".join(random.choices(string.ascii_letters, k=6)) for _ in range(int(1e4))]
text16 = np.random.choice(choices, n)
multitext = np.random.choice(["a", "b", "a; b"], n)
unique = [str(c) for c in range(n)]

df = pd.DataFrame(
    {
        "x": x,
        "y": y,
        "z": z,
        "text": text,
        "text16": text16,
        "multitext": multitext,
        "unique": unique,
    }
)

path = os.path.expanduser("~/mdv/dtypes_test")
p = MDVProject(path, delete_existing=True)
p.set_editable(True)
column_metadata = [
    {"name": "text", "datatype": "text"},
    {"name": "text16", "datatype": "text16"},
    {"name": "multitext", "datatype": "multitext"},
    {"name": "unique", "datatype": "unique"},
]
p.add_datasource("random data", df, columns=column_metadata)

# manually created & copied from debug chart dialog...
view = {
    "initialCharts": {
        "random data": [
            {
                "title": "x x y",
                "legend": "",
                "type": "wgl_scatter_plot",
                "param": ["x", "y"],
                "axis": {
                    "x": {
                        "size": 30,
                        "label": "x",
                        "textsize": 13,
                        "textSize": 13,
                        "tickfont": 10,
                    },
                    "y": {
                        "size": 45,
                        "label": "y",
                        "textsize": 13,
                        "textSize": 13,
                        "tickfont": 10,
                    },
                },
                "id": "g4lBjg",
                "size": [300, 300],
                "tooltip": {"show": False},
                "default_color": "#377eb8",
                "brush": "poly",
                "on_filter": "hide",
                "radius": 2.27,
                "opacity": 0.58,
                "color_by": "text",
                "color_legend": {"display": True, "pos": [45, 10]},
                "position": [0, 0],
            },
            {
                "title": "",
                "legend": "",
                "type": "selection_dialog_experimental",
                "param": ["multitext", "text", "text16", "unique", "x", "y", "z"],
                "filters": {
                    "multitext": {"category": ["a; b", "b", "a"]},
                    "text": None,
                    "text16": None,
                    "unique": None,
                    "x": [-4.647183, 5.0406637],
                    "y": None,
                    "z": None,
                },
                "id": "TEPaPn",
                "size": [539, 688],
                "position": [501, 0],
            },
            {
                "title": "",
                "legend": "",
                "type": "selection_dialog",
                "param": ["multitext", "text", "text16", "unique", "x", "y", "z"],
                "id": "k3y40k",
                "size": [583, 712],
                "filters": {
                    "multitext": {"category": [], "operand": "or"},
                    "text": {"category": [], "operand": None},
                    "x": None,
                    "y": None,
                    "z": None,
                },
                "position": [1059, 0],
            },
            {
                "title": "multitext",
                "legend": "",
                "type": "row_chart",
                "param": "multitext",
                "id": "ZBapEe",
                "size": [358, 508],
                "axis": {
                    "x": {"textSize": 13, "label": "", "size": 25, "tickfont": 10}
                },
                "position": [8, 298],
            },
        ]
    },
    "dataSources": {"random data": {"panelWidth": 100}},
}

p.set_view("default", view)
