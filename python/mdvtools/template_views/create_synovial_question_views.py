#!/usr/bin/env python3
"""Add synovial biological-question views without writing new columns.

Charts use columns and RNA rows that already exist in the project. Questions that
need a protein table, spatial coordinates, collection-score columns, or DEG
columns are left out, so datafile.h5 is not modified.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import h5py

from mdvtools.template_views.create_question_views import (
    HEAT_MAP_FEATURES,
    LinkedFeatures,
    atomic_write_json,
    box_plot,
    category_counts,
    cell_filter,
    composition_bar,
    compound_dot_plot,
    feature_list,
    feature_query,
    gene_info,
    heat_map,
    place_legends,
    placed,
    resolve_feature_fields,
    scatter,
    text_box,
)

DEFAULT_PROJECT = Path("~/mdv/mdv_synovial_latest")

CELLS = "cells"
RNA_DATASOURCE = "rna"
SUBGROUP = "rna_logged_counts"

UMAP_X = "rna:X_umap_mindist_0.25_1"
UMAP_Y = "rna:X_umap_mindist_0.25_2"
CELL_GROUP = "rna:bucket"
CELL_TYPE = "rna:annotations"
DIAGNOSIS = "rna:diagnosis"
SITE = "rna:site"
INFLAMMATION = "rna:inflammation_status"
PATIENT = "rna:patient_id"
HAD_TNF = "rna:had_tnf"
TNF_RESPONSE = "rna:tnf_response"
DATA_SET = "rna:data_set"

LAYOUT_WIDTH = 1590
DOT_PLOT_FEATURES = 5

MARKER_GENES = [
    "PTPRC", "CD3E", "CD4", "CXCL13", "PDCD1", "FOXP3", "KLRB1", "CD8A", "GZMK", "GZMB",
    "IFNG", "IL23R", "CD19", "MS4A1", "CD27", "MZB1", "DERL3", "IGHA1", "IGHG1", "CD14",
    "FCGR3A", "S100A8", "CD48", "C1QC", "IL1B", "MERTK", "LYVE1", "VSIG4", "TNF", "ETS2",
    "CD1C", "XCR1", "LAMP3", "FCER1A", "TCF4", "CPA3", "MKI67", "PDGFRA", "PRG4", "CLIC5",
    "THY1", "PDPN", "FAP", "POSTN", "COMP", "RSPO3", "C3", "SOX9", "CILP", "BGLAP", "CD34",
    "VWF", "PLVAP", "PODXL", "ACKR1", "TAGLN", "NOTCH3", "RGS5", "MYH11",
]

RNA = LinkedFeatures(datasource=RNA_DATASOURCE, subgroup=SUBGROUP, names=MARKER_GENES)

Q1_VIEW_NAME = "Where is my gene of interest expressed across cell populations?"
# Charts follow the genes picked in the genes panel (active links); Q1's panel starts on one gene.
Q1_START_GENE = "TNFRSF1B"
Q2_VIEW_NAME = "Is the expression different between disease stages/patient groups/tissue types?"
Q2_PATIENTS = ["CARTR02522", "CARTR02514", "CARTR02521", "CARTR02512", "CARTR02518"]
Q6_VIEW_NAME = "How comparable is the expression across different diseases?"
Q6_DIAGNOSES = ["RA", "PsA", "Healthy", "OA"]
# Q7 and Q9 were removed from the project. The names stay so older copies still sort and convert.
Q7_VIEW_NAME = "How comparable is the expression of gene collections/signatures across different diseases?"
Q8_VIEW_NAME = "Contrasting responders and non-responders to TNF inhibitor"
Q8_GROUPS = ["responder", "non_responder"]
Q8_GROUP_TITLES = {"responder": "Responder", "non_responder": "Non-responder"}
Q8_BOX_GENE = "TNF"
Q8_DATA_SETS = ["cartography"]
Q9_VIEW_NAME = "Comparing groups: is the UMAP different between TNF responders and non-responders?"

# Titles and axis labels only. Chart parameters keep the field id.
TITLE_BY_FIELD = {
    UMAP_X: "UMAP 1",
    UMAP_Y: "UMAP 2",
    "stroma:X_umap_1": "UMAP Stroma 1",
    "stroma:X_umap_2": "UMAP Stroma 2",
    "B_plasma:X_umap_1": "UMAP Bplasma 1",
    "B_plasma:X_umap_2": "UMAP Bplasma 2",
    "MNP:X_umap_1": "UMAP MNP 1",
    "MNP:X_umap_2": "UMAP MNP 2",
    "T_NK:X_umap_1": "UMAP T_NK 1",
    "T_NK:X_umap_2": "UMAP T_NK 2",
    "rna:annotations": "Cell Type",
    "rna:bucket": "Cell Group",
}
EXACT_TITLES = {
    "annotations": "Cell Type",
    "bucket": "Cell Group",
    "SJC28": "Swollen Joint Count (28)",
    "Sero_status": "Serostatus",
    "TJC28": "Tender Joint Count (28)",
    "sex": "Sex",
    "TNF Response": "TNF Treatment Response",
    "Had nbDMARD": "Prior nbDMARD",
    "Had TNF": "Prior TNF Inhibitor",
    "diagnosis": "Diagnosis",
    "sjc28": "Swollen Joint Count (28)",
    "esr": "ESR",
    "tjc28": "Tender Joint Count (28)",
    "disease_duration": "Disease Duration",
}
LEGACY_VIEWS = ("metadata", "Summary", "Stroma", "MNP", "T_NK", "B_plasma", "temp")


Q1_TEXT = f"""\
**Where is my gene of interest expressed across cell populations?** \
The first map shows the cell groups. The second map shows where the selected gene is expressed \
(darker = higher). In the **dot plot**, a bigger dot means more cells of that cell type express \
the gene, and a darker dot means a higher average level. The **heat map** shows the average level \
of the selected genes (up to {HEAT_MAP_FEATURES}) in each cell group. **MS4A1** in the gene list is CD20.

**To look at other genes:** click a gene in the list on the right (Ctrl/Cmd-click to pick \
several) and the charts update; the map shows the first gene picked.
"""

Q2_TEXT = f"""\
**Is the expression different between diseases, patients, tissues and inflammation states?** \
The **dot plots** compare the selected genes by **Diagnosis**, **Biopsy Site**, \
**Inflammation Status** and **Patient**: a bigger dot means more cells express the gene, and a \
darker dot means a higher average level.

**To compare groups:** pick one or more diagnoses, biopsy sites, cell types or inflammation states \
in **Choose groups to compare**, and every chart uses only those cells. **To change patients:** the **Patient** chart \
shows {len(Q2_PATIENTS)} example patients, a panel each, by cell group. Click its settings icon (⚙) \
and pick patients under **Split values**. **To change genes:** click genes in the list on the right \
(Ctrl/Cmd-click to pick several).
"""

Q6_TEXT = f"""\
**How comparable is the expression across different diseases?** The **heat map** shows the \
average level of the marker genes in every diagnosis, so genes that are high everywhere stand \
out from disease-specific ones. The **dot plot** shows the selected genes in each cell type, with \
a panel each for {", ".join(Q6_DIAGNOSES[:-1])} and {Q6_DIAGNOSES[-1]}: similar dots mean the gene \
behaves the same way in those diagnoses.

**To focus on a cell population:** pick cell types in **Focus on cell types** (all charts \
update). **To change genes in the dot plot:** click genes in the list on the right \
(Ctrl/Cmd-click to pick several). **To compare other diseases:** click the dot plot's settings icon \
(⚙) and choose diagnoses under **Split values**.
"""

Q8_TEXT = """\
**Contrasting TNF inhibitor responders and non-responders.** The **dot plot** shows the \
selected genes in each cell group, with a panel each for responders and non-responders. The **bar** shows the mix \
of cell groups, and the **box plot** shows how TNF varies in each response group. The page \
starts on RA samples from the cartography data set with a prior TNF inhibitor.

**To change the comparison:** choose a diagnosis, a prior TNF inhibitor value and a data set in \
**Choose patients to compare**. **To change genes:** click genes in the list on the right \
(Ctrl/Cmd-click to pick several).

**How the groups were defined:** Responder and Non-responder are the values of TNF Treatment \
Response. The page keeps cells whose Prior TNF Inhibitor value is Yes.
"""

def umap(title: str, color_by: str | dict[str, Any]) -> dict[str, Any]:
    """`color_by` is a field, or a `feature_query` for the selected gene (active link)."""
    return scatter(title, [UMAP_X, UMAP_Y], color_by, CELL_TYPE)


# Layouts are the final pixel positions from the hand-tuned synovial views (Nov 2026), so these
# views are not passed through spread_layout. Charts are [x, y], [width, height].


def marker_panel(start: list[str] | None = None) -> list[dict[str, Any]]:
    """Marker gene list (starting on `start`, default all markers) above the Gene Info chart."""
    return [
        placed(feature_list("Marker genes", start or MARKER_GENES), [0, 0], [270, 400]),
        placed(gene_info("Click a gene to plot it"), [0, 425], [270, 455]),
    ]


def rna_view(
    cells: list[dict[str, Any]], panel: list[dict[str, Any]] | None = None
) -> dict[str, Any]:
    return {
        "dataSources": {
            CELLS: {"panelWidth": 85},
            RNA_DATASOURCE: {"panelWidth": 15},
        },
        "initialCharts": {CELLS: cells, RNA_DATASOURCE: panel or marker_panel()},
    }


def build_q1_view(fields: dict[str, str]) -> dict[str, Any]:
    del fields
    return rna_view(
        [
            placed(text_box("About this page", Q1_TEXT), [0, 0], [LAYOUT_WIDTH, 120]),
            placed(umap("Cell Group", CELL_GROUP), [84, 143], [585, 340]),
            placed(umap("Selected Gene Expression", feature_query(RNA, 1)), [794, 142], [587, 344]),
            placed(
                compound_dot_plot("Selected genes by Cell Type", CELL_TYPE, feature_query(RNA, DOT_PLOT_FEATURES)),
                [4, 533],
                [929, 733],
            ),
            placed(
                heat_map("Cell Group", CELL_GROUP, feature_query(RNA, HEAT_MAP_FEATURES), y_label_size=108),
                [981, 522],
                [626, 628],
            ),
        ],
        marker_panel([Q1_START_GENE]),
    )


def build_q2_view(fields: dict[str, str]) -> dict[str, Any]:
    del fields
    selected = feature_query(RNA, DOT_PLOT_FEATURES)
    return rna_view([
        placed(text_box("About this page", Q2_TEXT), [0, 0], [1100, 200]),
        placed(compound_dot_plot("Selected genes by Diagnosis", DIAGNOSIS, selected), [0, 225], [482, 655]),
        # Dot plots are sized to leave room for their legends above the plot.
        placed(compound_dot_plot("Selected genes by Biopsy Site", SITE, selected), [506, 225], [521, 655]),
        placed(
            compound_dot_plot("Selected genes by Inflammation Status", INFLAMMATION, selected),
            [1050, 230],
            [558, 415],
        ),
        placed(
            cell_filter("Choose groups to compare", [DIAGNOSIS, SITE, CELL_TYPE, INFLAMMATION]),
            [1051, 697],
            [517, 511],
        ),
        # The React dot plot can't limit its rows to a few patients, so it has a panel per
        # example patient instead.
        placed(
            compound_dot_plot(
                "Selected genes by Patient",
                CELL_GROUP,
                selected,
                split=(PATIENT, Q2_PATIENTS),
                y_label_size=110,
            ),
            [0, 905],
            [1013, 385],
        ),
    ])


def build_q6_view(fields: dict[str, str]) -> dict[str, Any]:
    return rna_view([
        placed(text_box("About this page", Q6_TEXT), [0, 0], [1547, 150]),
        placed(
            heat_map(
                "Marker genes by Diagnosis (average expression)",
                DIAGNOSIS,
                [fields[gene] for gene in MARKER_GENES],
                y_label_size=168,
                y_tickfont=7,
            ),
            [0, 175],
            [771, 809],
        ),
        placed(
            compound_dot_plot(
                "Selected genes by Diagnosis",
                CELL_TYPE,
                feature_query(RNA, DOT_PLOT_FEATURES),
                split=(DIAGNOSIS, Q6_DIAGNOSES),
            ),
            [787, 175],
            [795, 907],
        ),
        placed(cell_filter("Focus on cell types", [CELL_TYPE]), [30, 1001], [734, 299]),
    ])


def build_q8_view(fields: dict[str, str]) -> dict[str, Any]:
    return rna_view([
        placed(text_box("About this page", Q8_TEXT), [0, 0], [1093, 191]),
        placed(
            cell_filter(
                "Choose patients to compare",
                [DIAGNOSIS, HAD_TNF, TNF_RESPONSE, DATA_SET],
                {DIAGNOSIS: ["RA"], HAD_TNF: ["Yes"], TNF_RESPONSE: Q8_GROUPS, DATA_SET: Q8_DATA_SETS},
            ),
            [1117, 0],
            [473, 190],
        ),
        placed(
            compound_dot_plot(
                "Selected genes by TNF Treatment Response",
                CELL_GROUP,
                feature_query(RNA, DOT_PLOT_FEATURES),
                split=(TNF_RESPONSE, Q8_GROUPS),
                y_label_size=110,
            ),
            [0, 220],
            [878, 485],
        ),
        placed(
            box_plot(f"{Q8_BOX_GENE} by TNF Treatment Response", TNF_RESPONSE, fields[Q8_BOX_GENE]),
            [916, 212],
            [672, 296],
        ),
        placed(
            composition_bar("Cell Group mix by TNF Treatment Response", TNF_RESPONSE, CELL_GROUP),
            [922, 541],
            [662, 310],
        ),
    ])


VIEW_BUILDERS = {
    Q1_VIEW_NAME: build_q1_view,
    Q2_VIEW_NAME: build_q2_view,
    Q6_VIEW_NAME: build_q6_view,
    Q8_VIEW_NAME: build_q8_view,
}


def datasource(datasources: list[dict[str, Any]], name: str) -> dict[str, Any]:
    found = next((item for item in datasources if item.get("name") == name), None)
    if found is None:
        raise ValueError(f"Project has no '{name}' datasource.")
    return found


def validate_project(datasources: list[dict[str, Any]]) -> None:
    cells = datasource(datasources, CELLS)
    fields = {column.get("field") for column in cells.get("columns", [])}
    required = [
        UMAP_X, UMAP_Y, CELL_GROUP, CELL_TYPE, DIAGNOSIS, SITE, INFLAMMATION, PATIENT,
        HAD_TNF, TNF_RESPONSE, DATA_SET,
    ]
    missing = [field for field in required if field not in fields]
    if missing:
        raise ValueError(f"'{CELLS}' is missing fields: {missing}")
    link = cells.get("links", {}).get(RNA_DATASOURCE, {})
    subgroups = link.get("rows_as_columns", {}).get("subgroups", {})
    if SUBGROUP not in subgroups:
        raise ValueError(f"No '{SUBGROUP}' rows-as-columns link from '{CELLS}' to '{RNA_DATASOURCE}'.")
    gene_names = set(
        next(column["values"] for column in datasource(datasources, RNA_DATASOURCE)["columns"] if column.get("field") == "name")
    )
    absent = [gene for gene in [*RNA.names, Q1_START_GENE] if gene not in gene_names]
    if absent:
        raise ValueError(f"Genes missing from '{RNA_DATASOURCE}': {absent}")
    for field, categories in (
        (DIAGNOSIS, Q6_DIAGNOSES),
        (DATA_SET, Q8_DATA_SETS),
        (HAD_TNF, ["Yes"]),
        (TNF_RESPONSE, Q8_GROUPS),
        (PATIENT, Q2_PATIENTS),
    ):
        values = next(column["values"] for column in cells["columns"] if column.get("field") == field)
        unknown = [category for category in categories if category not in values]
        if unknown:
            raise ValueError(f"'{field}' has no categories {unknown}")


def retitle_text(value: str) -> str:
    exact = EXACT_TITLES.get(value)
    if exact is None:
        exact = EXACT_TITLES.get(value.strip())
    if exact is not None:
        return exact
    updated = value
    for field, title in sorted(TITLE_BY_FIELD.items(), key=lambda item: len(item[0]), reverse=True):
        updated = updated.replace(field, title)
    return updated


def use_tnk_umap(view: dict[str, Any]) -> int:
    """The T_NK maps were stored on the MNP UMAP. Point them at the T_NK UMAP."""
    changed = 0
    charts = (view.get("initialCharts") or {}).get(CELLS) or []
    if not isinstance(charts, list):
        return 0
    for chart in charts:
        if not isinstance(chart, dict) or chart.get("param") != ["MNP:X_umap_1", "MNP:X_umap_2"]:
            continue
        chart["param"] = ["T_NK:X_umap_1", "T_NK:X_umap_2"]
        axis = chart.get("axis")
        if not isinstance(axis, dict):
            axis = {}
            chart["axis"] = axis
        for key, label in (("x", "UMAP T_NK 1"), ("y", "UMAP T_NK 2")):
            spec = axis.get(key)
            if not isinstance(spec, dict):
                axis[key] = {"label": label}
            else:
                spec["label"] = label
        changed += 1
    return changed


def retitle_existing_views(views: dict[str, Any], only: set[str] | None = None) -> int:
    """Point stored titles at the sheet's final names. Field ids in param stay put."""
    changed = 0
    for name, view in views.items():
        if only is not None and name not in only:
            continue
        if not isinstance(view, dict):
            continue
        if name == "T_NK":
            changed += use_tnk_umap(view)
        for charts in (view.get("initialCharts") or {}).values():
            if not isinstance(charts, list):
                continue
            for chart in charts:
                if not isinstance(chart, dict):
                    continue
                title = chart.get("title")
                if isinstance(title, str):
                    updated = retitle_text(title)
                    if updated != title:
                        chart["title"] = updated
                        changed += 1
                axis = chart.get("axis")
                if isinstance(axis, dict):
                    for axis_spec in axis.values():
                        if isinstance(axis_spec, dict) and isinstance(axis_spec.get("label"), str):
                            updated = retitle_text(axis_spec["label"])
                            if updated != axis_spec["label"]:
                                axis_spec["label"] = updated
                                changed += 1
    return changed


def prepend_text(view: dict[str, Any], markdown: str) -> None:
    """Put `markdown` above the text of the view's first text box."""
    for charts in view["initialCharts"].values():
        for chart in charts:
            if chart.get("type") == "text_box_chart":
                chart["text"] = f"{markdown.strip()}\n\n{chart['text']}"
                return


def h5_stamp(project_path: Path) -> tuple[int, int]:
    stat = (project_path / "datafile.h5").stat()
    return stat.st_size, stat.st_mtime_ns


def add_synovial_question_views(
    project_path: Path, *, dry_run: bool = False, treatment_disclaimer: str | None = None
) -> dict[str, Any]:
    """`treatment_disclaimer` (Markdown) goes above the text of the TNF response page."""
    project_path = project_path.expanduser().resolve()
    views_path = project_path / "views.json"
    state_path = project_path / "state.json"
    datasources_path = project_path / "datasources.json"
    for path in (views_path, state_path, datasources_path, project_path / "datafile.h5"):
        if not path.is_file():
            raise FileNotFoundError(f"Not a complete MDV project; missing {path}")

    before = h5_stamp(project_path)
    datasources = json.loads(datasources_path.read_text(encoding="utf-8"))
    validate_project(datasources)
    with h5py.File(project_path / "datafile.h5", "r") as handle:
        fields = resolve_feature_fields(handle, datasources, RNA)
    after_read = h5_stamp(project_path)
    if after_read != before:
        raise RuntimeError(f"Reading datafile.h5 changed it: {before} -> {after_read}")

    views = json.loads(views_path.read_text(encoding="utf-8"))
    state = json.loads(state_path.read_text(encoding="utf-8"))
    retitled = retitle_existing_views(views)
    all_views = list(state.get("all_views") or [])
    counts = category_counts(datasources)
    for name, build in VIEW_BUILDERS.items():
        views[name] = place_legends(build(fields), counts)
        if name not in all_views:
            all_views.append(name)
    state["all_views"] = all_views
    if treatment_disclaimer:
        prepend_text(views[Q8_VIEW_NAME], treatment_disclaimer)

    summary: dict[str, Any] = {
        "project": str(project_path),
        "views": list(VIEW_BUILDERS),
        "skipped": [
            "protein concordance (no protein datasource)",
            "spatial localisation (no spatial coordinates)",
            "collection score columns",
            "DEG columns",
        ],
        "retitled_charts": retitled,
        "h5_bytes": before[0],
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
    after_write = h5_stamp(project_path)
    if after_write != before:
        raise RuntimeError(f"datafile.h5 changed while writing views: {before} -> {after_write}")
    summary["backups"] = backups
    return summary


def update_legacy_views(project_path: Path, *, dry_run: bool = False) -> dict[str, Any]:
    """Retitle the non-question views and point T_NK at its own UMAP. Question views stay as stored."""
    project_path = project_path.expanduser().resolve()
    views_path = project_path / "views.json"
    state_path = project_path / "state.json"
    if not views_path.is_file() or not state_path.is_file():
        raise FileNotFoundError(f"Not a complete MDV project; missing views or state in {project_path}")
    h5_path = project_path / "datafile.h5"
    before = h5_stamp(project_path) if h5_path.is_file() else None
    views = json.loads(views_path.read_text(encoding="utf-8"))
    state = json.loads(state_path.read_text(encoding="utf-8"))
    question_names = set(VIEW_BUILDERS)
    before_questions = {name: json.dumps(views[name], sort_keys=True) for name in question_names if name in views}
    changed = retitle_existing_views(views, only=set(LEGACY_VIEWS))
    after_questions = {name: json.dumps(views[name], sort_keys=True) for name in question_names if name in views}
    if before_questions != after_questions:
        raise RuntimeError("Legacy update changed a question view.")
    missing = [name for name in LEGACY_VIEWS if name not in views]
    if missing:
        raise ValueError(f"Legacy views missing from the project: {missing}")
    all_views = list(state.get("all_views") or [])
    for name in list(views):
        if name not in all_views:
            all_views.append(name)
    summary: dict[str, Any] = {
        "project": str(project_path),
        "legacy_views": list(LEGACY_VIEWS),
        "question_views": [name for name in all_views if name in question_names],
        "changed": changed,
        "dry_run": dry_run,
    }
    if dry_run:
        return summary
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    backup_path = views_path.with_name(f"{views_path.name}.bak.{timestamp}")
    backup_path.write_bytes(views_path.read_bytes())
    atomic_write_json(views_path, views)
    if before is not None and h5_stamp(project_path) != before:
        raise RuntimeError("datafile.h5 changed while updating legacy views.")
    summary["backup"] = str(backup_path)
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", type=Path, default=DEFAULT_PROJECT)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument(
        "--treatment-disclaimer",
        type=Path,
        help="Markdown file to show above the text of the TNF response page.",
    )
    parser.add_argument(
        "--legacy-only",
        action="store_true",
        help="Retitle the existing non-question views and fix the T_NK UMAP. Do not rewrite question views.",
    )
    args = parser.parse_args()
    if args.legacy_only:
        print(json.dumps(update_legacy_views(args.project, dry_run=args.dry_run), indent=2))
        return 0
    disclaimer = args.treatment_disclaimer.expanduser().read_text(encoding="utf-8") if args.treatment_disclaimer else None
    summary = add_synovial_question_views(args.project, dry_run=args.dry_run, treatment_disclaimer=disclaimer)
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
