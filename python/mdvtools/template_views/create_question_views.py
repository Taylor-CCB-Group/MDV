#!/usr/bin/env python3
"""Add biological-question views to an MDV cell atlas project.

Each view answers one question for a non-specialist user: a short markdown guide,
marker gene / protein lists on the genes / proteins panels, and cell charts that follow
the selected features through the cells -> genes and cells -> proteins rows-as-columns
links. Layouts are absolute so a whole view fits on one screen.
"""

from __future__ import annotations

import argparse
import json
import os
import random
import string
import tempfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import h5py
import numpy as np
from scipy import stats

from mdvtools.mdvproject import MDVProject

DEFAULT_PROJECT = Path("~/mdv/myeloid_v2_testing_fields_mutated")

CELLS = "cells"
NAME_FIELD = "name"


@dataclass(frozen=True)
class LinkedFeatures:
    datasource: str
    subgroup: str
    names: list[str]


GENES = LinkedFeatures(
    datasource="genes",
    subgroup="gs",
    names=[
        "CD14", "S100A8", "TNFRSF14", "ITGAL", "CD274", "ICAM1", "FCGR3A", "CD48", "CD68", "MERTK",
        "LYVE1", "ELMO1", "TNF", "IL1B", "MX1", "CD1C", "TCF4", "XCR1", "LAMP3", "FCER1A",
        "ETS2", "SPP1", "GPNMB", "TREM2", "TREM1", "VSIG4", "IL3RA", "CPA3", "CD34", "MKI67",
    ],
)
PROTEINS = LinkedFeatures(
    datasource="proteins",
    subgroup="pro_rc",
    names=[
        "adt_CD48", "adt_CD14", "adt_CLEC12A", "adt_CD11b", "adt_CD11a", "adt_CD16", "adt_CX3CR1",
        "adt_CD1c", "adt_TNFRSF14", "adt_PD-L1", "adt_CD83", "adt_ICAM1", "adt_CD40", "adt_CD86",
        "adt_HLA-DR", "adt_CD33", "adt_CD13", "adt_CD141", "adt_IL3RA", "adt_CR1", "adt_CD38",
        "adt_CD82",
    ],
)

UMAP_X = "X_umap_mindist_0.25_1"
UMAP_Y = "X_umap_mindist_0.25_2"
PROTEIN_UMAP_X = "prot:X_umap_1"
CELL_COMPARTMENT = "bucket_caf"
CELL_SUBGROUP = "sub_bucket_caf"
CELL_STATE = "annotation_caf"

DIAGNOSIS = "diagnosis"
DISEASE = "disease"
PATIENT = "Patient"
TISSUE = "tissue_simple"
INFLAMMATION = "inflammation_status"
BATCH = "batch"

DOT_PLOT_FEATURES = 5

CELLS_PANEL_WIDTH = 85
GENES_PANEL_WIDTH = 15
WIDE_LAYOUT_WIDTH = 1590

PROTEIN_VIEW_PANEL_WIDTHS = {CELLS: 70, GENES.datasource: 15, PROTEINS.datasource: 15}
PROTEIN_VIEW_LAYOUT_WIDTH = 1310

Q1_VIEW_NAME = "Where is my gene of interest expressed across cell populations?"
Q1_FEATURE_GENES = GENES.names[:3]

Q1_TEXT = f"""\
**Where is my gene of interest expressed across cell populations?** \
The first map shows the myeloid cell subgroups. The next maps show where \
{", ".join(Q1_FEATURE_GENES)} are expressed (darker = higher). In the **dot plot**, a bigger dot \
means more cells in that cell state express the gene, and a darker dot means a higher average \
level. The **heat map** shows the average level of all {len(GENES.names)} marker genes in each cell state.

**To look at other genes:** click a gene in the list on the right (Ctrl/Cmd-click to pick \
several) and the dot plot updates. Type a gene name in **Marker genes** to add it to the list.
"""


Q2_VIEW_NAME = "Is the expression different between disease stages/patient groups/tissue types?"
Q2_PATIENTS = ["CARTG00030", "CARTD01030", "CARTR02006", "CARTG00006", "CARTR02014"]

Q2_TEXT = f"""\
**Is the expression different between diseases, patients, tissues and inflammation states?** \
The **dot plots** compare the selected genes by **Diagnosis**, **Disease**, **Tissue**, \
**Inflammation Status** and **Patient**: a bigger dot means more cells express the gene, and a \
darker dot means a higher average level.

**To compare groups:** pick one or more diagnoses, diseases or cell states in **Choose groups to \
compare**, and every chart uses only those cells. **To change patients:** the **Patient** chart \
shows {len(Q2_PATIENTS)} example patients. Click its settings icon (⚙) and pick patients under \
**Background Filter Categories**. **To change genes:** click genes in the list on the right \
(Ctrl/Cmd-click to pick several).
"""


Q3_VIEW_NAME = "Is the expression concordant with protein expression (CITE-seq)?"
Q3_PAIRS = [("CD14", "adt_CD14"), ("FCGR3A", "adt_CD16")]
Q3_MATCHING = [
    ("CD14", "adt_CD14"), ("FCGR3A", "adt_CD16"), ("ITGAL", "adt_CD11a"), ("CD274", "adt_PD-L1"),
    ("ICAM1", "adt_ICAM1"), ("TNFRSF14", "adt_TNFRSF14"), ("CD48", "adt_CD48"),
    ("CD1C", "adt_CD1c"), ("IL3RA", "adt_IL3RA"),
]

Q3_TEXT = f"""\
**Is the expression concordant with protein expression?** Protein was measured (CITE-seq) for \
most cells. The left dot plot shows **RNA** for the selected genes and the right one shows \
**protein** for the selected proteins, in each cell state. Both use only batches with protein \
data. The maps below compare {" and ".join(f"{g} RNA with {p}" for g, p in Q3_PAIRS)} protein: \
similar colour patterns mean RNA and protein agree.

**To compare a gene with its protein:** click the gene in the genes list and its protein in the \
proteins list (Ctrl/Cmd-click to pick several). Matching pairs: \
{", ".join(f"{g} = {p}" for g, p in Q3_MATCHING)}. Protein values are raw counts, so compare \
colour (level) rather than dot size.
"""


TOTAL_COUNTS = "total_counts"
GENE_COLLECTION_FIELD = "gene_collection"
NO_GENE_COLLECTION = "Other"
GENE_COLLECTIONS = {
    "Monocyte": ["CD14", "S100A8", "FCGR3A", "ITGAL", "CD48", "TNFRSF14"],
    "Tissue macrophage": ["CD68", "MERTK", "LYVE1", "VSIG4", "ELMO1"],
    "Lipid-associated macrophage": ["SPP1", "GPNMB", "TREM2"],
    "Inflammatory response": ["TNF", "IL1B", "TREM1", "ICAM1", "CD274", "MX1", "ETS2"],
    "Dendritic cell": ["CD1C", "FCER1A", "XCR1", "LAMP3", "TCF4", "IL3RA"],
    "Progenitor and proliferating": ["CD34", "CPA3", "MKI67"],
}


def collection_score_field(collection: str) -> str:
    return "score_" + collection.lower().replace(" ", "_").replace("-", "_")


Q6_VIEW_NAME = "How comparable is the expression across different diseases?"
Q6_DIAGNOSES = ["Crohns disease", "Ulcerative colitis"]
Q6_BOX_GENE = GENES.names[0]

Q6_TEXT = f"""\
**How comparable is the expression across different diseases?** The **heat map** shows the \
average level of all marker genes in every diagnosis, so genes that are high everywhere stand \
out from disease-specific ones. The **box plot** shows how {Q6_BOX_GENE} varies between cells \
within each diagnosis. The **dot plot** shows the selected genes in each cell state, with a \
panel each for {" and ".join(Q6_DIAGNOSES)}: similar dots mean the gene behaves the same way in \
both diseases.

**To focus on a cell population:** pick cell states in **Focus on cell states** (all charts \
update). **To change genes in the dot plots:** click genes in the list on the right \
(Ctrl/Cmd-click to pick several). **To compare other diseases:** click the dot plot's settings icon (⚙) \
and choose diagnoses under **Split values**.
"""

Q7_VIEW_NAME = "How comparable is the expression of gene collections/signatures across different diseases?"
Q7_DEFAULT_COLLECTION = next(iter(GENE_COLLECTIONS))
Q7_DIAGNOSES = Q6_DIAGNOSES

Q7_TEXT = f"""\
**How comparable are gene collections across different diseases?** A gene collection \
(signature) is a group of marker genes that together mark a cell type or process. Its **score** \
is the average normalised expression of its genes in each cell. The **heat map** shows the average \
score of every collection in each diagnosis, so collections that are high in only some diseases \
stand out. The **box plot** shows how the {Q7_DEFAULT_COLLECTION} score varies between cells within \
each diagnosis. The **dot plot** shows the individual genes of the chosen collection in each \
cell state, with a panel each for {" and ".join(Q7_DIAGNOSES)}.

**To choose a collection:** pick one in **Choose a gene collection** on the right; the dot plots \
then show its genes. **To focus on a cell population:** pick cell states in **Focus on cell \
states** (all charts update). **To compare other diseases:** click the dot plot's settings icon \
(⚙) and choose diagnoses under **Split values**.
"""


TREATMENT = "treatment"
RESPONSE = "response"
RESPONSE_CLASS = "response_class"
TREATMENT_STAGE = "treatment_stage"
THERAPY = "therapy"
RESPONSE_STAGE = "response_stage"
NOT_DETERMINED = "ND"
RESPONDER = "R"
NON_RESPONDER = "NR"
RESPONSE_CLASSES = {
    "responsive": RESPONDER,
    "good-response": RESPONDER,
    "good_response": RESPONDER,
    "moderate-response": RESPONDER,
    "Remission": RESPONDER,
    "non-responsive": NON_RESPONDER,
    "poor-response": NON_RESPONDER,
    "minor-response": NON_RESPONDER,
    "Non_Remission": NON_RESPONDER,
}
TREATMENT_STAGES = {"Pre": "baseline", "Post": "post-treatment"}
THERAPIES = {"TNF": "anti-TNF", "DUPI": "Dupilumab", "Uste": "Ustekinumab", "Vedo": "Vedolizumab"}
RESPONSE_STAGE_GROUPS = [
    f"{response} {stage}"
    for stage in TREATMENT_STAGES.values()
    for response in (RESPONDER, NON_RESPONDER)
]

Q8_VIEW_NAME = "Contrasting R vs NR at baseline and at post treatment"
Q8_DEFAULTS = {DIAGNOSIS: ["Crohns disease"], TISSUE: ["Gut"], THERAPY: ["anti-TNF"]}
Q8_BOX_COLLECTION = "Inflammatory response"

Q8_TEXT = f"""\
**Contrasting responders (R) and non-responders (NR) at baseline and post-treatment.** The \
**dot plot** shows the selected genes in each cell subgroup, with one panel per group: R and NR at \
baseline, then R and NR post-treatment. Compare neighbouring panels for R vs NR, and the first \
pair with the second for the effect of treatment. The **heat map** shows the gene collection \
scores in the four groups, and the **box plot** shows the {Q8_BOX_COLLECTION} score per group. The page starts on \
Crohn's disease gut biopsies treated with anti-TNF.

**To change the comparison:** choose a disease, tissue and therapy in **Choose patients to \
compare** (pick one of each so like is compared with like). **To change genes:** click genes in \
the list on the right (Ctrl/Cmd-click to pick several).

**How R and NR were defined:** baseline = "Pre_" and post-treatment = "Post_" samples in \
Treatment. R = responsive, good response, moderate response (EULAR convention) or remission. \
NR = non-responsive, poor response, minor response or non-remission. Cells with an unknown \
response are left out.
"""

DEG_COMPARISON = Q8_DEFAULTS
DEG_COMPARISON_LABEL = "Crohn's disease gut biopsies on anti-TNF"
DEG_GROUPS = ("baseline", "post-treatment")
DEG_MIN_FRACTION = 0.1
DEG_MIN_LOG2_FOLD = 0.5
DEG_MAX_P = 0.05
DEG_TOP = 5
DEG_LOG2_FOLD = "deg_log2_fold_change"
DEG_EVIDENCE = "deg_neg_log10_p"
DEG_ADJUSTED_P = "deg_adjusted_p"
DEG_CHANGE = "deg_change"
DEG_TOP_FIELD = "deg_top"
DEG_UP = "Up after treatment (significant)"
DEG_DOWN = "Down after treatment (significant)"
DEG_CANDIDATE_UP = "Up after treatment (candidate)"
DEG_CANDIDATE_DOWN = "Down after treatment (candidate)"
DEG_NO_CHANGE = "No clear change"
DEG_NOT_TESTED = "Not tested (rarely expressed)"
DEG_TOP_UP = f"Top {DEG_TOP} up"
DEG_TOP_DOWN = f"Top {DEG_TOP} down"
DEG_COLUMN_NAMES = {
    DEG_LOG2_FOLD: "Change after treatment (log2 fold)",
    DEG_EVIDENCE: "Evidence (-log10 p)",
    DEG_ADJUSTED_P: "Adjusted p",
    DEG_CHANGE: "Change after treatment",
    DEG_TOP_FIELD: "Top changed genes",
}

Q9_VIEW_NAME = "Comparing groups: is the UMAP different between treatment states? DEG?"
Q9_PANEL_WIDTHS = {CELLS: 75, GENES.datasource: 25}

Q9_TEXT = f"""\
**Is the UMAP different between treatment states, and which genes change?** The two **UMAPs** \
show the same cells before (left) and after (right) treatment, coloured by cell subgroup; a \
region that fills or empties points to a shift in cell types. The **bar** shows the proportion of \
each cell subgroup at each stage. The **dot plot** shows the selected genes in each cell subgroup, \
with a panel for baseline and one for post-treatment.

**Differentially expressed genes (DEG)** on the right compare each patient's post-treatment \
sample with their own baseline sample, in {DEG_COMPARISON_LABEL} (the starting selection). In \
the **volcano plot**, genes to the right go up after treatment, genes to the left go down, and \
higher means stronger evidence. The page starts on the top {DEG_TOP} candidates in each \
direction. **To explore other genes:** clear **Show genes**, then draw round genes in the volcano \
plot or click them in the table (Ctrl/Cmd-click to pick several).

**Good to know:** in this comparison no gene changes clearly once all genes tested are taken \
into account, so the "candidate" genes (p < {DEG_MAX_P} and at least {2 ** DEG_MIN_LOG2_FOLD:.1f}-fold, \
before correction) are leads to follow up rather than findings; "significant" would mean adjusted \
p < {DEG_MAX_P}. Changing **Choose patients to compare** updates the charts on the left, but the \
DEG list stays fixed to {DEG_COMPARISON_LABEL}. A change in the cell mix can also make a gene look \
changed.
"""


def chart_id() -> str:
    return "".join(random.choices(string.ascii_letters + string.digits, k=6))


def feature_field(features: LinkedFeatures, name: str, row_index: int) -> str:
    return f"{features.subgroup}|{name} ({features.subgroup})|{row_index}"


def feature_query(features: LinkedFeatures, max_items: int) -> dict[str, Any]:
    return {
        "type": "RowsAsColsQuery",
        "linkedDsName": features.datasource,
        "maxItems": max_items,
        "subgroupName": features.subgroup,
    }


def placed(chart: dict[str, Any], position: list[int], size: list[int]) -> dict[str, Any]:
    return {**chart, "position": position, "size": size, "id": chart_id()}


def text_box(title: str, text: str) -> dict[str, Any]:
    return {"type": "text_box_chart", "title": title, "param": [], "text": text}


def umap(
    title: str, color_by: str, only: tuple[str, str] | None = None
) -> dict[str, Any]:
    """`only=(column, value)` restricts this chart (not the others) to cells with that value."""
    chart = scatter(title, [UMAP_X, UMAP_Y], color_by, CELL_STATE)
    if only:
        chart["category_filters"] = [{"column": only[0], "category": [only[1]]}]
    return chart


def scatter(title: str, xy: list[str], color_by: str, tooltip: str) -> dict[str, Any]:
    return {
        "type": "wgl_scatter_plot_dev",
        "title": title,
        "legend": "",
        "param": xy,
        "color_by": color_by,
        "color_legend": {"display": True, "pos": [5, 5]},
        "field_legend": {"display": True},
        "tooltip": {"show": True, "column": [tooltip]},
        "dimension": "2d",
        "point_shape": "circle",
        "radius": 10,
        "course_radius": 1,
        "opacity": 1,
        "on_filter": "grey",
        "zoom_on_filter": False,
        "axis": {
            "x": {"rotate_labels": False, "size": 30, "tickfont": 9},
            "y": {"rotate_labels": False, "size": 35, "tickfont": 9},
        },
    }


def dot_plot(
    title: str,
    category: str,
    features: str | list[str] | dict[str, Any],
    y_label_size: int = 175,
    only: tuple[str, list[str]] | None = None,
) -> dict[str, Any]:
    """`only=(column, values)` restricts this chart (not the others) to cells with those values."""
    feature_params = features if isinstance(features, list) else [features]
    return {
        "type": "dot_plot",
        "title": title,
        "legend": "",
        "param": [category, *feature_params],
        "background_filter": {
            "column": only[0] if only else None,
            "categories": only[1] if only else [],
        },
        "axis": {
            "x": {"label": "", "rotate_labels": True, "size": 70, "textSize": 13, "tickfont": 10},
            "y": {"label": "", "size": y_label_size, "textSize": 13, "tickfont": 10},
            "ry": {"label": "", "size": 25, "textSize": 13, "tickfont": 10},
        },
        "color_scale": {"log": False, "palette": "reds"},
        "color_legend": {"display": True, "pos": [40, 10]},
        "fraction_legend": {"display": True, "pos": [0, 0], "position": [0, 0]},
        "y_axis_order": "data",
    }


def compound_dot_plot(
    title: str,
    category: str,
    features: str | list[str] | dict[str, Any],
    split: tuple[str, list[str]],
    y_label_size: int = 175,
) -> dict[str, Any]:
    """One dot plot with a panel per value in `split=(column, values)`, all on one colour and size scale.

    Use it in place of several `dot_plot`s that differ only in the value they are restricted to.
    """
    feature_params = features if isinstance(features, list) else [features]
    return {
        "type": "compound_dot_plot",
        "title": title,
        "legend": "",
        "param": [category, *feature_params],
        "split_by": split[0],
        "split_values": list(split[1]),
        "axis": {
            "x": {
                "label": "",
                "rotate_labels": True,
                "size": 70,
                "textSize": 13,
                "tickfont": 10,
                "fit_labels": True,
                "truncate_labels": True,
                "thin_labels": False,
            },
            "y": {"label": "", "rotate_labels": False, "size": y_label_size, "textSize": 13, "tickfont": 10},
        },
        "color_scale": {"log": False, "palette": "reds", "mode": "auto"},
        "color_legend": {"display": True},
        "fraction_legend": {"display": True},
        "reserve_legend_space": True,
        "y_axis_order": "data",
    }


def heat_map(
    title: str, category: str, fields: list[str], y_label_size: int = 90
) -> dict[str, Any]:
    return {
        "type": "heat_map",
        "title": title,
        "legend": "",
        "param": [category, *fields],
        "method": "mean",
        "cluster_rows": False,
        "cluster_columns": False,
        "axis": {
            "x": {"label": "", "rotate_labels": True, "size": 150, "textSize": 13, "tickfont": 9},
            "y": {"label": "", "size": y_label_size, "textSize": 13, "tickfont": 9},
        },
        "color_scale": {"log": False, "trim": "0.05"},
        "color_legend": {"display": True, "pos": [0, 0]},
    }


def box_plot(title: str, category: str, field: str) -> dict[str, Any]:
    return {
        "type": "box_plot",
        "title": title,
        "legend": "",
        "param": [category, field],
        "color_by": category,
        "color_legend": {"display": False},
        "axis": {
            "x": {"label": "", "rotate_labels": True, "size": 110, "textSize": 13, "tickfont": 10},
            "y": {"label": "Expression", "size": 45, "textSize": 13, "tickfont": 10},
        },
        "tooltip": {"show": False},
        "default_color": "#377eb8",
        "brush": "poly",
        "on_filter": "hide",
        "radius": 2,
        "opacity": 0.5,
        "band_width": 0,
        "intervals": 20,
    }


def composition_bar(title: str, group_by: str, split_by: str) -> dict[str, Any]:
    return {
        "type": "stacked_row_chart",
        "title": title,
        "legend": "",
        "param": [group_by, split_by],
        "sort": "name",
        "color_legend": {"display": True, "pos": [5, 5]},
        "axis": {
            "x": {"label": "Proportion of cells", "size": 35, "textSize": 13, "tickfont": 10},
            "y": {"label": "", "size": 110, "textSize": 13, "tickfont": 11},
        },
    }


def category_list(title: str, field: str, selected: list[str]) -> dict[str, Any]:
    return {
        "type": "selection_dialog",
        "title": title,
        "legend": "",
        "param": [field],
        "filters": {field: {"category": selected}},
        "order": {},
        "noClearFilters": False,
    }


def cell_filter(
    title: str, fields: list[str], selected: dict[str, list[str]] | None = None
) -> dict[str, Any]:
    selected = selected or {}
    return {
        "type": "selection_dialog",
        "title": title,
        "legend": "",
        "param": fields,
        "filters": {
            field: {"category": selected[field]} if field in selected else None for field in fields
        },
        "order": {},
        "noClearFilters": False,
    }


def feature_list(title: str, names: list[str]) -> dict[str, Any]:
    return {
        "type": "selection_dialog",
        "title": title,
        "legend": "",
        "param": [NAME_FIELD],
        "filters": {NAME_FIELD: {"category": names}},
        "order": {},
        "noClearFilters": False,
    }


def feature_table(title: str, extra_columns: tuple[str, ...] = ()) -> dict[str, Any]:
    return {
        "type": "table_chart",
        "title": title,
        "legend": "",
        "param": [NAME_FIELD, *extra_columns],
        "column_widths": {},
        "include_index": False,
    }


def genes_panel() -> list[dict[str, Any]]:
    return [
        placed(feature_list("Marker genes", GENES.names), [0, 0], [270, 400]),
        placed(feature_table("Click a gene to plot it"), [0, 405], [270, 455]),
    ]


def proteins_panel() -> list[dict[str, Any]]:
    return [
        placed(feature_list("Marker proteins", PROTEINS.names), [0, 0], [270, 400]),
        placed(feature_table("Click a protein to plot it"), [0, 405], [270, 455]),
    ]


def collections_panel() -> list[dict[str, Any]]:
    return [
        placed(
            category_list("Choose a gene collection", GENE_COLLECTION_FIELD, [Q7_DEFAULT_COLLECTION]),
            [0, 0],
            [270, 300],
        ),
        placed(
            feature_table("Genes in the collection", (GENE_COLLECTION_FIELD,)), [0, 305], [270, 555]
        ),
    ]


def genes_view(
    cells: list[dict[str, Any]], genes: list[dict[str, Any]] | None = None
) -> dict[str, Any]:
    return {
        "dataSources": {
            CELLS: {"panelWidth": CELLS_PANEL_WIDTH},
            GENES.datasource: {"panelWidth": GENES_PANEL_WIDTH},
        },
        "initialCharts": {CELLS: cells, GENES.datasource: genes or genes_panel()},
    }


@dataclass(frozen=True)
class ProjectFields:
    features: dict[str, str]
    protein_batches: list[str]


def build_q1_view(fields: ProjectFields) -> dict[str, Any]:
    features = [
        placed(umap(gene, fields.features[gene]), [435 + i * 385, 125], [380, 350])
        for i, gene in enumerate(Q1_FEATURE_GENES)
    ]
    return genes_view([
        placed(text_box("About this page", Q1_TEXT), [0, 0], [WIDE_LAYOUT_WIDTH, 120]),
        placed(umap("Cell Subgroup", CELL_SUBGROUP), [0, 125], [430, 350]),
        *features,
        placed(
            dot_plot(
                "Selected genes by Cell State", CELL_STATE, feature_query(GENES, DOT_PLOT_FEATURES)
            ),
            [0, 480],
            [600, 380],
        ),
        placed(
            heat_map(
                "Marker genes by Cell State (average expression)",
                CELL_STATE,
                [fields.features[g] for g in GENES.names],
            ),
            [605, 480],
            [985, 380],
        ),
    ])


def build_q2_view(fields: ProjectFields) -> dict[str, Any]:
    selected = feature_query(GENES, DOT_PLOT_FEATURES)
    return genes_view([
        placed(text_box("About this page", Q2_TEXT), [0, 0], [1100, 200]),
        placed(
            cell_filter("Choose groups to compare", [DIAGNOSIS, DISEASE, CELL_STATE]),
            [1105, 0],
            [485, 200],
        ),
        placed(dot_plot("Selected genes by Diagnosis", DIAGNOSIS, selected), [0, 205], [470, 655]),
        placed(dot_plot("Selected genes by Disease", DISEASE, selected), [475, 205], [555, 325]),
        placed(dot_plot("Selected genes by Tissue", TISSUE, selected), [475, 535], [555, 325]),
        placed(
            dot_plot("Selected genes by Inflammation Status", INFLAMMATION, selected),
            [1035, 205],
            [555, 325],
        ),
        placed(
            dot_plot(
                "Selected genes by Patient",
                PATIENT,
                selected,
                y_label_size=90,
                only=(PATIENT, Q2_PATIENTS),
            ),
            [1035, 535],
            [555, 325],
        ),
    ])


def build_q3_view(fields: ProjectFields) -> dict[str, Any]:
    with_protein = (BATCH, fields.protein_batches)
    half = PROTEIN_VIEW_LAYOUT_WIDTH // 2
    quarter = PROTEIN_VIEW_LAYOUT_WIDTH // 4
    pair_maps = []
    for i, (gene, protein) in enumerate(Q3_PAIRS):
        for j, (name, kind) in enumerate(((gene, "RNA"), (protein, "protein"))):
            pair_maps.append(
                placed(
                    umap(f"{name} ({kind})", fields.features[name]),
                    [(2 * i + j) * quarter, 505],
                    [quarter - 5, 355],
                )
            )
    cells = [
        placed(text_box("About this page", Q3_TEXT), [0, 0], [PROTEIN_VIEW_LAYOUT_WIDTH, 150]),
        placed(
            dot_plot(
                "RNA: selected genes by Cell State",
                CELL_STATE,
                feature_query(GENES, DOT_PLOT_FEATURES),
                only=with_protein,
            ),
            [0, 155],
            [half - 5, 345],
        ),
        placed(
            dot_plot(
                "Protein: selected proteins by Cell State",
                CELL_STATE,
                feature_query(PROTEINS, DOT_PLOT_FEATURES),
                only=with_protein,
            ),
            [half, 155],
            [half, 345],
        ),
        *pair_maps,
    ]
    return {
        "dataSources": {ds: {"panelWidth": width} for ds, width in PROTEIN_VIEW_PANEL_WIDTHS.items()},
        "initialCharts": {
            CELLS: cells,
            GENES.datasource: genes_panel(),
            PROTEINS.datasource: proteins_panel(),
        },
    }


def build_q6_view(fields: ProjectFields) -> dict[str, Any]:
    selected = feature_query(GENES, DOT_PLOT_FEATURES)
    per_disease = placed(
        compound_dot_plot(
            f"Selected genes in {' and '.join(Q6_DIAGNOSES)}",
            CELL_STATE,
            selected,
            split=(DIAGNOSIS, Q6_DIAGNOSES),
        ),
        [795, 505],
        [795, 355],
    )
    return genes_view([
        placed(text_box("About this page", Q6_TEXT), [0, 0], [1100, 150]),
        placed(cell_filter("Focus on cell states", [CELL_STATE]), [1105, 0], [485, 150]),
        placed(
            heat_map(
                "Marker genes by Diagnosis (average expression)",
                DIAGNOSIS,
                [fields.features[g] for g in GENES.names],
            ),
            [0, 155],
            [790, 705],
        ),
        placed(
            box_plot(f"{Q6_BOX_GENE} by Diagnosis", DIAGNOSIS, fields.features[Q6_BOX_GENE]),
            [795, 155],
            [795, 345],
        ),
        per_disease,
    ])


def build_q7_view(fields: ProjectFields) -> dict[str, Any]:
    collection_genes = feature_query(GENES, max(len(g) for g in GENE_COLLECTIONS.values()))
    per_disease = placed(
        compound_dot_plot(
            f"Collection genes in {' and '.join(Q7_DIAGNOSES)}",
            CELL_STATE,
            collection_genes,
            split=(DIAGNOSIS, Q7_DIAGNOSES),
        ),
        [795, 155],
        [795, 705],
    )
    return genes_view(
        [
            placed(text_box("About this page", Q7_TEXT), [0, 0], [1100, 150]),
            placed(cell_filter("Focus on cell states", [CELL_STATE]), [1105, 0], [485, 150]),
            placed(
                heat_map(
                    "Gene collection scores by Diagnosis (average)",
                    DIAGNOSIS,
                    [collection_score_field(c) for c in GENE_COLLECTIONS],
                    y_label_size=190,
                ),
                [0, 155],
                [790, 400],
            ),
            placed(
                box_plot(
                    f"{Q7_DEFAULT_COLLECTION} score by Diagnosis",
                    DIAGNOSIS,
                    collection_score_field(Q7_DEFAULT_COLLECTION),
                ),
                [0, 560],
                [790, 300],
            ),
            per_disease,
        ],
        collections_panel(),
    )


def build_q8_view(fields: ProjectFields) -> dict[str, Any]:
    selected = feature_query(GENES, DOT_PLOT_FEATURES)
    groups = placed(
        compound_dot_plot(
            "Selected genes by group",
            CELL_SUBGROUP,
            selected,
            split=(RESPONSE_STAGE, RESPONSE_STAGE_GROUPS),
            y_label_size=90,
        ),
        [0, 220],
        [1185, 640],
    )
    return genes_view([
        placed(text_box("About this page", Q8_TEXT), [0, 0], [1100, 215]),
        placed(
            cell_filter(
                "Choose patients to compare",
                [DIAGNOSIS, TISSUE, THERAPY, RESPONSE_STAGE],
                {**Q8_DEFAULTS, RESPONSE_STAGE: RESPONSE_STAGE_GROUPS},
            ),
            [1105, 0],
            [485, 215],
        ),
        groups,
        placed(
            heat_map(
                "Gene collection scores by group",
                RESPONSE_STAGE,
                [collection_score_field(c) for c in GENE_COLLECTIONS],
                y_label_size=190,
            ),
            [1190, 220],
            [400, 318],
        ),
        placed(
            box_plot(
                f"{Q8_BOX_COLLECTION} score by group",
                RESPONSE_STAGE,
                collection_score_field(Q8_BOX_COLLECTION),
            ),
            [1190, 543],
            [400, 317],
        ),
    ])


def deg_panel() -> list[dict[str, Any]]:
    return [
        placed(
            category_list("Show genes", DEG_TOP_FIELD, [DEG_TOP_UP, DEG_TOP_DOWN]), [0, 0], [465, 150]
        ),
        placed(
            scatter(
                "Volcano plot: change after treatment", [DEG_LOG2_FOLD, DEG_EVIDENCE], DEG_CHANGE, NAME_FIELD
            ),
            [0, 155],
            [465, 400],
        ),
        placed(
            feature_table("Genes (click to plot)", (DEG_CHANGE, DEG_LOG2_FOLD, DEG_ADJUSTED_P)),
            [0, 560],
            [465, 300],
        ),
    ]


def build_q9_view(fields: ProjectFields) -> dict[str, Any]:
    selected = feature_query(GENES, 2 * DEG_TOP)
    umaps = [
        placed(
            umap(f"Cell Subgroup at {stage}", CELL_SUBGROUP, only=(TREATMENT_STAGE, stage)),
            [i * 465, 205],
            [460, 400],
        )
        for i, stage in enumerate(DEG_GROUPS)
    ]
    by_stage = placed(
        compound_dot_plot(
            f"Selected genes at {' and '.join(DEG_GROUPS)}",
            CELL_SUBGROUP,
            selected,
            split=(TREATMENT_STAGE, list(DEG_GROUPS)),
            y_label_size=90,
        ),
        [0, 610],
        [1400, 250],
    )
    return {
        "dataSources": {ds: {"panelWidth": width} for ds, width in Q9_PANEL_WIDTHS.items()},
        "initialCharts": {
            CELLS: [
                placed(text_box("About this page", Q9_TEXT), [0, 0], [900, 200]),
                placed(
                    cell_filter(
                        "Choose patients to compare",
                        [DIAGNOSIS, TISSUE, THERAPY, TREATMENT_STAGE],
                        {**DEG_COMPARISON, TREATMENT_STAGE: list(DEG_GROUPS)},
                    ),
                    [905, 0],
                    [495, 200],
                ),
                *umaps,
                placed(
                    composition_bar("Cell mix by treatment stage", TREATMENT_STAGE, CELL_SUBGROUP),
                    [930, 205],
                    [470, 400],
                ),
                by_stage,
            ],
            GENES.datasource: deg_panel(),
        },
    }


VIEW_BUILDERS = {
    Q1_VIEW_NAME: build_q1_view,
    Q2_VIEW_NAME: build_q2_view,
    Q3_VIEW_NAME: build_q3_view,
    Q6_VIEW_NAME: build_q6_view,
    Q7_VIEW_NAME: build_q7_view,
    Q8_VIEW_NAME: build_q8_view,
    Q9_VIEW_NAME: build_q9_view,
}


def datasource(datasources: list[dict[str, Any]], name: str) -> dict[str, Any]:
    found = next((ds for ds in datasources if ds.get("name") == name), None)
    if found is None:
        raise ValueError(f"Project has no '{name}' datasource.")
    return found


def text_values(datasources: list[dict[str, Any]], ds_name: str, field: str) -> list[str]:
    column = next(c for c in datasource(datasources, ds_name)["columns"] if c.get("field") == field)
    return column["values"]


def resolve_feature_fields(
    handle: h5py.File, datasources: list[dict[str, Any]], features: LinkedFeatures
) -> dict[str, str]:
    """Map feature names to rows-as-columns field names using their h5 row index.

    The index in the field name is the row in the linked datasource, which is not the
    position of the name in the column's `values` list.
    """
    row_of = feature_rows(handle, datasources, features)
    return {n: feature_field(features, n, row_of[n]) for n in features.names}


def row_names(handle: h5py.File, datasources: list[dict[str, Any]], ds_name: str) -> list[str]:
    values = text_values(datasources, ds_name, NAME_FIELD)
    return [values[code] for code in handle[ds_name][NAME_FIELD][:]]


def feature_rows(
    handle: h5py.File, datasources: list[dict[str, Any]], features: LinkedFeatures
) -> dict[str, int]:
    row_of: dict[str, int] = {}
    for row, name in enumerate(row_names(handle, datasources, features.datasource)):
        row_of.setdefault(name, row)
    missing = [n for n in features.names if n not in row_of]
    if missing:
        raise ValueError(f"Not found in '{features.datasource}': {missing}")
    return row_of


def compute_collection_scores(
    handle: h5py.File, datasources: list[dict[str, Any]]
) -> dict[str, np.ndarray]:
    """Per cell, the mean over a collection's genes of log1p(counts per 10k)."""
    cells = datasource(datasources, CELLS)
    subgroup = cells["links"][GENES.datasource]["rows_as_columns"]["subgroups"][GENES.subgroup]
    matrix = handle[CELLS][subgroup["name"]]
    pointers = matrix["p"][:]
    rows = feature_rows(handle, datasources, GENES)
    scale = 1e4 / handle[CELLS][TOTAL_COUNTS][:]
    scores: dict[str, np.ndarray] = {}
    for collection, genes in GENE_COLLECTIONS.items():
        total = np.zeros(cells["size"], dtype=np.float64)
        for gene in genes:
            start, end = pointers[rows[gene]], pointers[rows[gene] + 1]
            expressing = matrix["i"][start:end]
            total[expressing] += np.log1p(matrix["x"][start:end] * scale[expressing])
        scores[collection] = (total / len(genes)).astype(np.float32)
    return scores


def response_columns(responses: list[str], treatments: list[str]) -> dict[str, list[str]]:
    """Harmonised R/NR, treatment stage (Pre_/Post_ prefix), therapy and their combination."""
    classes = [RESPONSE_CLASSES.get(r, NOT_DETERMINED) for r in responses]
    stages: list[str] = []
    therapies: list[str] = []
    for treatment in treatments:
        prefix, _, drug = treatment.partition("_")
        known = prefix in TREATMENT_STAGES and drug in THERAPIES
        stages.append(TREATMENT_STAGES[prefix] if known else NOT_DETERMINED)
        therapies.append(THERAPIES[drug] if known else NOT_DETERMINED)
    combined = [
        f"{c} {s}" if NOT_DETERMINED not in (c, s) else NOT_DETERMINED
        for c, s in zip(classes, stages)
    ]
    return {RESPONSE_CLASS: classes, TREATMENT_STAGE: stages, THERAPY: therapies, RESPONSE_STAGE: combined}


RESPONSE_COLUMN_NAMES = {
    RESPONSE_CLASS: "Responder (R/NR)",
    TREATMENT_STAGE: "Treatment Stage",
    THERAPY: "Therapy",
    RESPONSE_STAGE: "Response and Stage",
}


def has_columns(
    project_path: Path, datasources: list[dict[str, Any]], wanted: list[tuple[str, str]]
) -> bool:
    """Whether every (datasource, field) is both described in datasources.json and stored in the h5.

    A shared datasources.json can describe columns that a copy's older datafile.h5 lacks.
    """
    described = {(ds["name"], c["field"]) for ds in datasources for c in ds.get("columns", [])}
    with h5py.File(project_path / "datafile.h5", "r") as handle:
        return all(key in described and key[1] in handle[key[0]] for key in wanted)


def add_response_columns(project_path: Path, datasources: list[dict[str, Any]]) -> list[str]:
    """Add the derived response columns to cells; returns new fields."""
    if has_columns(project_path, datasources, [(CELLS, f) for f in RESPONSE_COLUMN_NAMES]):
        return []
    with h5py.File(project_path / "datafile.h5", "r") as handle:
        codes = {f: handle[CELLS][f][:] for f in (RESPONSE, TREATMENT)}
    labels = {f: text_values(datasources, CELLS, f) for f in (RESPONSE, TREATMENT)}
    columns = response_columns(
        [labels[RESPONSE][c] for c in codes[RESPONSE]],
        [labels[TREATMENT][c] for c in codes[TREATMENT]],
    )
    backup_datasources(project_path)
    project = MDVProject(str(project_path))
    for field, values in columns.items():
        project.set_column(
            CELLS, {"name": RESPONSE_COLUMN_NAMES[field], "field": field, "datatype": "text"}, values
        )
    return list(columns)


def cell_labels(handle: h5py.File, datasources: list[dict[str, Any]], field: str) -> np.ndarray:
    return np.array(text_values(datasources, CELLS, field), dtype=object)[handle[CELLS][field][:]]


def compute_treatment_deg(
    handle: h5py.File, datasources: list[dict[str, Any]]
) -> dict[str, np.ndarray]:
    """Paired pseudobulk DEG, post-treatment vs baseline.

    Cells are averaged on log1p(counts per 10k) per patient and stage, and the per-patient
    post - baseline differences are tested (one-sample t-test), so patients with both samples
    are the replicates.
    """
    cells = datasource(datasources, CELLS)
    in_comparison = np.ones(cells["size"], dtype=bool)
    for field, values in DEG_COMPARISON.items():
        in_comparison &= np.isin(cell_labels(handle, datasources, field), values)
    stage = cell_labels(handle, datasources, TREATMENT_STAGE)
    patient = handle[CELLS][PATIENT][:].astype(np.int64)
    in_stage = [in_comparison & (stage == name) for name in DEG_GROUPS]
    paired = np.intersect1d(*(np.unique(patient[cells]) for cells in in_stage))
    sample = np.full(cells["size"], -1, dtype=np.int64)
    for k, code in enumerate(paired):
        for g, cells_in_stage in enumerate(in_stage):
            sample[cells_in_stage & (patient == code)] = 2 * k + g
    n_samples = 2 * len(paired)
    sample_sizes = np.bincount(sample[sample >= 0], minlength=n_samples).astype(np.float64)
    groups = np.arange(n_samples) % 2

    subgroup = cells["links"][GENES.datasource]["rows_as_columns"]["subgroups"][GENES.subgroup]
    matrix = handle[CELLS][subgroup["name"]]
    pointers = matrix["p"][:]
    n_genes = len(pointers) - 1
    scale = 1e4 / handle[CELLS][TOTAL_COUNTS][:]
    sample_sums = np.zeros((n_genes, n_samples))
    sample_expressing = np.zeros((n_genes, n_samples))
    block = 500
    for first in range(0, n_genes, block):
        last = min(first + block, n_genes)
        start, end = pointers[first], pointers[last]
        cell_index = matrix["i"][start:end]
        gene = np.repeat(np.arange(first, last), np.diff(pointers[first : last + 1]))
        s = sample[cell_index]
        keep = s >= 0
        value = np.log1p(matrix["x"][start:end][keep] * scale[cell_index[keep]])
        bins = (gene[keep] - first) * n_samples + s[keep]
        size = (last - first) * n_samples
        sample_sums[first:last] = np.bincount(bins, value, size).reshape(-1, n_samples)
        sample_expressing[first:last] = np.bincount(bins, minlength=size).reshape(-1, n_samples)

    sample_means = sample_sums / sample_sizes
    group_cells = np.array([sample_sizes[groups == g].sum() for g in range(2)])
    fraction = np.stack(
        [sample_expressing[:, groups == g].sum(axis=1) / group_cells[g] for g in range(2)], axis=1
    )
    baseline, post = sample_means[:, groups == 0], sample_means[:, groups == 1]
    difference = post - baseline
    tested = fraction.max(axis=1) >= DEG_MIN_FRACTION
    tested &= difference.std(axis=1) > 0
    p = np.full(n_genes, np.nan)
    p[tested] = stats.ttest_1samp(difference[tested], 0, axis=1).pvalue
    adjusted = np.full(n_genes, np.nan)
    adjusted[tested] = stats.false_discovery_control(p[tested])
    fold = np.log2(np.expm1(post.mean(axis=1)) + 1e-9) - np.log2(np.expm1(baseline.mean(axis=1)) + 1e-9)

    change = np.full(n_genes, DEG_NOT_TESTED, dtype=object)
    change[tested] = DEG_NO_CHANGE
    large_up = tested & (fold >= DEG_MIN_LOG2_FOLD)
    large_down = tested & (fold <= -DEG_MIN_LOG2_FOLD)
    nominal = np.nan_to_num(p, nan=1.0) < DEG_MAX_P
    significant = np.nan_to_num(adjusted, nan=1.0) < DEG_MAX_P
    change[large_up & nominal] = DEG_CANDIDATE_UP
    change[large_down & nominal] = DEG_CANDIDATE_DOWN
    change[large_up & significant] = DEG_UP
    change[large_down & significant] = DEG_DOWN
    top = np.full(n_genes, "Other", dtype=object)
    for label, direction in ((DEG_TOP_UP, large_up), (DEG_TOP_DOWN, large_down)):
        candidates = np.flatnonzero(direction & nominal)
        top[candidates[np.argsort(p[candidates])][:DEG_TOP]] = label
    return {
        DEG_LOG2_FOLD: np.where(tested, fold, np.nan).astype(np.float32),
        DEG_EVIDENCE: (-np.log10(p)).astype(np.float32),
        DEG_ADJUSTED_P: adjusted.astype(np.float32),
        DEG_CHANGE: change,
        DEG_TOP_FIELD: top,
    }


def add_treatment_deg(project_path: Path, datasources: list[dict[str, Any]]) -> list[str]:
    """Add DEG result columns to genes; returns new fields."""
    if has_columns(project_path, datasources, [(GENES.datasource, f) for f in DEG_COLUMN_NAMES]):
        return []
    with h5py.File(project_path / "datafile.h5", "r") as handle:
        columns = compute_treatment_deg(handle, datasources)
    backup_datasources(project_path)
    project = MDVProject(str(project_path))
    for field, values in columns.items():
        datatype = "text" if values.dtype == object else "double"
        project.set_column(
            GENES.datasource,
            {"name": DEG_COLUMN_NAMES[field], "field": field, "datatype": datatype},
            list(values) if datatype == "text" else values,
        )
    return list(columns)


def backup_datasources(project_path: Path) -> None:
    datasources_path = project_path / "datasources.json"
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    datasources_path.with_name(f"{datasources_path.name}.bak.{timestamp}").write_bytes(
        datasources_path.read_bytes()
    )


def add_gene_collections(project_path: Path, datasources: list[dict[str, Any]]) -> list[str]:
    """Add collection score columns to cells and a collection column to genes; returns new fields."""
    score_fields = [collection_score_field(c) for c in GENE_COLLECTIONS]
    wanted = [(CELLS, f) for f in score_fields] + [(GENES.datasource, GENE_COLLECTION_FIELD)]
    if has_columns(project_path, datasources, wanted):
        return []

    with h5py.File(project_path / "datafile.h5", "r") as handle:
        scores = compute_collection_scores(handle, datasources)
        gene_names = row_names(handle, datasources, GENES.datasource)
    collection_of = {g: c for c, genes in GENE_COLLECTIONS.items() for g in genes}

    backup_datasources(project_path)
    project = MDVProject(str(project_path))
    for collection, score in scores.items():
        project.set_column(
            CELLS,
            {"name": f"{collection} score", "field": collection_score_field(collection), "datatype": "double"},
            score,
        )
    project.set_column(
        GENES.datasource,
        {"name": "Gene Collection", "field": GENE_COLLECTION_FIELD, "datatype": "text"},
        [collection_of.get(name, NO_GENE_COLLECTION) for name in gene_names],
    )
    return [f for _, f in wanted]


def resolve_protein_batches(handle: h5py.File, datasources: list[dict[str, Any]]) -> list[str]:
    """Batches containing any cell with protein data (a protein UMAP coordinate)."""
    values = text_values(datasources, CELLS, BATCH)
    has_protein = np.isfinite(handle[CELLS][PROTEIN_UMAP_X][:])
    codes = np.unique(handle[CELLS][BATCH][:][has_protein])
    return [values[c] for c in codes]


def resolve_project_fields(project_path: Path, datasources: list[dict[str, Any]]) -> ProjectFields:
    with h5py.File(project_path / "datafile.h5", "r") as handle:
        features = {
            **resolve_feature_fields(handle, datasources, GENES),
            **resolve_feature_fields(handle, datasources, PROTEINS),
        }
        return ProjectFields(features, resolve_protein_batches(handle, datasources))


def validate_fields(datasources: list[dict[str, Any]]) -> None:
    cells = datasource(datasources, CELLS)
    cell_fields = {column.get("field") for column in cells.get("columns", [])}
    required = (
        UMAP_X, UMAP_Y, PROTEIN_UMAP_X, CELL_COMPARTMENT, CELL_SUBGROUP, CELL_STATE,
        DIAGNOSIS, DISEASE, PATIENT, TISSUE, INFLAMMATION, BATCH, TOTAL_COUNTS, RESPONSE, TREATMENT,
    )
    outside_markers = {g for genes in GENE_COLLECTIONS.values() for g in genes} - set(GENES.names)
    if outside_markers:
        raise ValueError(f"Gene collections use genes outside the marker list: {outside_markers}")
    missing = [f for f in required if f not in cell_fields]
    if missing:
        raise ValueError(f"'{CELLS}' is missing fields: {missing}")
    for features in (GENES, PROTEINS):
        datasource(datasources, features.datasource)
        link = cells.get("links", {}).get(features.datasource, {})
        if features.subgroup not in link.get("rows_as_columns", {}).get("subgroups", {}):
            raise ValueError(
                f"No '{features.subgroup}' rows-as-columns link from '{CELLS}' to '{features.datasource}'."
            )


def atomic_write_json(path: Path, data: Any) -> None:
    fd, temporary_path = tempfile.mkstemp(prefix=f"{path.name}.", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(data, handle, indent=2)
            handle.write("\n")
        os.replace(temporary_path, path)
    except Exception:
        if os.path.exists(temporary_path):
            os.unlink(temporary_path)
        raise


def add_question_views(project_path: Path, *, dry_run: bool = False) -> dict[str, Any]:
    project_path = project_path.expanduser().resolve()
    views_path = project_path / "views.json"
    state_path = project_path / "state.json"
    datasources_path = project_path / "datasources.json"
    for path in (views_path, state_path, datasources_path):
        if not path.is_file():
            raise FileNotFoundError(f"Not a complete MDV project; missing {path}")

    datasources = json.loads(datasources_path.read_text(encoding="utf-8"))
    validate_fields(datasources)
    added_fields: list[str] = []
    if not dry_run:
        added_fields = add_gene_collections(project_path, datasources)
        datasources = json.loads(datasources_path.read_text(encoding="utf-8"))
        added_fields += add_response_columns(project_path, datasources)
        datasources = json.loads(datasources_path.read_text(encoding="utf-8"))
        added_fields += add_treatment_deg(project_path, datasources)
        datasources = json.loads(datasources_path.read_text(encoding="utf-8"))
    fields = resolve_project_fields(project_path, datasources)
    views = json.loads(views_path.read_text(encoding="utf-8"))
    state = json.loads(state_path.read_text(encoding="utf-8"))

    all_views = list(state.get("all_views") or [])
    for name, build in VIEW_BUILDERS.items():
        views[name] = build(fields)
        if name not in all_views:
            all_views.append(name)
    state["all_views"] = all_views

    summary: dict[str, Any] = {
        "project": str(project_path),
        "views": list(VIEW_BUILDERS),
        "protein_batches": fields.protein_batches,
        "added_fields": added_fields,
        "dry_run": dry_run,
    }
    if dry_run:
        return summary

    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    backups = []
    for path in (views_path, state_path):
        backup_path = path.with_name(f"{path.name}.bak.{timestamp}")
        backup_path.write_bytes(path.read_bytes())
        backups.append(str(backup_path))
    atomic_write_json(views_path, views)
    atomic_write_json(state_path, state)
    summary["backups"] = backups
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", type=Path, default=DEFAULT_PROJECT, help="MDV project directory.")
    parser.add_argument("--dry-run", action="store_true", help="Validate without writing files.")
    args = parser.parse_args()
    print(json.dumps(add_question_views(args.project, dry_run=args.dry_run), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
