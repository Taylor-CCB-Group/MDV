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
    LinkedFeatures,
    atomic_write_json,
    box_plot,
    cell_filter,
    composition_bar,
    dot_plot,
    feature_list,
    feature_query,
    feature_table,
    heat_map,
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

# Same collections as the myeloid atlas, without IL3RA, which is not in this RNA table.
GENE_COLLECTIONS = {
    "Monocyte": ["CD14", "S100A8", "FCGR3A", "ITGAL", "CD48", "TNFRSF14"],
    "Tissue macrophage": ["CD68", "MERTK", "LYVE1", "VSIG4", "ELMO1"],
    "Lipid-associated macrophage": ["SPP1", "GPNMB", "TREM2"],
    "Inflammatory response": ["TNF", "IL1B", "TREM1", "ICAM1", "CD274", "MX1", "ETS2"],
    "Dendritic cell": ["CD1C", "FCER1A", "XCR1", "LAMP3", "TCF4"],
    "Progenitor and proliferating": ["CD34", "CPA3", "MKI67"],
}

RNA = LinkedFeatures(datasource=RNA_DATASOURCE, subgroup=SUBGROUP, names=MARKER_GENES)
COLLECTION_GENES = list(dict.fromkeys(gene for genes in GENE_COLLECTIONS.values() for gene in genes))
PLOTTED_GENES = LinkedFeatures(
    datasource=RNA_DATASOURCE,
    subgroup=SUBGROUP,
    names=list(dict.fromkeys([*MARKER_GENES, *COLLECTION_GENES])),
)

Q1_VIEW_NAME = "Where is my gene of interest expressed across cell populations?"
Q1_FEATURE_GENES = ["CD14", "S100A8", "TNF"]
Q2_VIEW_NAME = "Is the expression different between disease stages/patient groups/tissue types?"
Q2_PATIENTS = ["CARTR02522", "CARTR02514", "CARTR02521", "CARTR02512", "CARTR02518"]
Q6_VIEW_NAME = "How comparable is the expression across different diseases?"
Q6_DIAGNOSES = ["RA", "PsA"]
Q6_BOX_GENE = "CD14"
Q7_VIEW_NAME = "How comparable is the expression of gene collections/signatures across different diseases?"
Q8_VIEW_NAME = "Contrasting responders and non-responders to TNF inhibitor"
Q8_GROUPS = ["responder", "non_responder"]
Q8_GROUP_TITLES = {"responder": "Responder", "non_responder": "Non-responder"}
Q8_BOX_GENE = "TNF"
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
}


def collection_sentence() -> str:
    names = list(GENE_COLLECTIONS)
    return ", ".join(names[:-1]) + f", and {names[-1]}"


Q1_TEXT = f"""\
**Where is my gene of interest expressed across cell populations?** \
The first map shows the cell groups. The next maps show where \
{", ".join(Q1_FEATURE_GENES)} are expressed (darker = higher). In the **dot plot**, a bigger dot \
means more cells of that cell type express the gene, and a darker dot means a higher average \
level. The **heat map** shows the average level of the marker genes in each cell group. \
**MS4A1** in the gene list is CD20.

**To look at other genes:** click a gene in the list on the right (Ctrl/Cmd-click to pick \
several) and the dot plot updates.
"""

Q2_TEXT = f"""\
**Is the expression different between diseases, patients, tissues and inflammation states?** \
The **dot plots** compare the selected genes by **Diagnosis**, **Biopsy Site**, \
**Inflammation Status** and **Patient**: a bigger dot means more cells express the gene, and a \
darker dot means a higher average level.

**To compare groups:** pick one or more diagnoses, biopsy sites or cell types in **Choose groups \
to compare**, and every chart uses only those cells. **To change patients:** the **Patient** chart \
shows {len(Q2_PATIENTS)} example patients. Click its settings icon (⚙) and pick patients under \
**Background Filter Categories**. **To change genes:** click genes in the list on the right \
(Ctrl/Cmd-click to pick several).
"""

Q6_TEXT = f"""\
**How comparable is the expression across different diseases?** The **heat map** shows the \
average level of the marker genes in every diagnosis, so genes that are high everywhere stand \
out from disease-specific ones. The **box plot** shows how {Q6_BOX_GENE} varies between cells \
within each diagnosis. The two **dot plots** show the selected genes in each cell type for \
{" and ".join(Q6_DIAGNOSES)}: similar dots mean the gene behaves the same way in both diseases.

**To focus on a cell population:** pick cell types in **Focus on cell types** (all charts \
update). **To change genes in the dot plots:** click genes in the list on the right \
(Ctrl/Cmd-click to pick several). **To change a dot plot's disease:** click its settings icon (⚙) \
and choose a diagnosis under **Background Filter Categories**.
"""

Q7_TEXT = f"""\
**How comparable are gene collections across different diseases?** A gene collection is a group \
of marker genes that together mark a cell type or process. The collections are {collection_sentence()}. \
The **heat map** shows the average level of each of those genes in every diagnosis. The two \
**dot plots** show the genes you select, in each cell type, for {" and ".join(Q6_DIAGNOSES)}.

**To choose genes:** click genes in the list on the right (Ctrl/Cmd-click to pick several). \
**To focus on a cell population:** pick cell types in **Focus on cell types**. **To change a dot \
plot's disease:** click its settings icon (⚙) and choose a diagnosis under **Background Filter \
Categories**.
"""

Q8_TEXT = """\
**Contrasting TNF inhibitor responders and non-responders.** The two **dot plots** show the \
selected genes in each cell group for responders and non-responders. The **bar** shows the mix \
of cell groups, and the **box plot** shows how TNF varies in each response group. The page \
starts on RA samples with a prior TNF inhibitor.

**To change the comparison:** choose a diagnosis and a prior TNF inhibitor value in **Choose \
patients to compare**. **To change genes:** click genes in the list on the right \
(Ctrl/Cmd-click to pick several).

**How the groups were defined:** Responder and Non-responder are the values of TNF Treatment \
Response. The page keeps cells whose Prior TNF Inhibitor value is Yes.
"""

Q9_TEXT = """\
**Is the UMAP different between TNF responders and non-responders?** The two **UMAPs** show \
responder cells and non-responder cells, coloured by cell group. A region that fills or empties \
points to a shift in cell groups. The **bar** shows the proportion of each cell group in the two \
response groups. The **dot plots** show the selected genes in each cell group. The page starts \
on RA samples with a prior TNF inhibitor.

**To change the comparison:** choose a diagnosis and a prior TNF inhibitor value in **Choose \
patients to compare**. **To change genes:** click genes in the list on the right \
(Ctrl/Cmd-click to pick several). Response is TNF Treatment Response for cells with Prior TNF \
Inhibitor set to Yes.
"""


def umap(title: str, color_by: str, only: tuple[str, str] | None = None) -> dict[str, Any]:
    chart = scatter(title, [UMAP_X, UMAP_Y], color_by, CELL_TYPE)
    if only:
        chart["category_filters"] = [{"column": only[0], "category": [only[1]]}]
    return chart


def marker_panel(names: list[str]) -> list[dict[str, Any]]:
    return [
        placed(feature_list("Marker genes", names), [0, 0], [270, 400]),
        placed(feature_table("Click a gene to plot it"), [0, 405], [270, 455]),
    ]


def rna_view(
    cells: list[dict[str, Any]], panel: list[dict[str, Any]] | None = None
) -> dict[str, Any]:
    return {
        "dataSources": {
            CELLS: {"panelWidth": 85},
            RNA_DATASOURCE: {"panelWidth": 15},
        },
        "initialCharts": {CELLS: cells, RNA_DATASOURCE: panel or marker_panel(MARKER_GENES)},
    }


def build_q1_view(fields: dict[str, str]) -> dict[str, Any]:
    maps = [
        placed(umap(gene, fields[gene]), [390 + index * 400, 125], [395, 340])
        for index, gene in enumerate(Q1_FEATURE_GENES)
    ]
    return rna_view([
        placed(text_box("About this page", Q1_TEXT), [0, 0], [LAYOUT_WIDTH, 120]),
        placed(umap("Cell Group", CELL_GROUP), [0, 125], [385, 340]),
        *maps,
        placed(
            dot_plot("Selected genes by Cell Type", CELL_TYPE, feature_query(RNA, DOT_PLOT_FEATURES)),
            [0, 475],
            [700, 390],
        ),
        placed(
            heat_map(
                "Marker genes by Cell Group (average expression)",
                CELL_GROUP,
                [fields[gene] for gene in MARKER_GENES],
            ),
            [710, 475],
            [880, 390],
        ),
    ])


def build_q2_view(fields: dict[str, str]) -> dict[str, Any]:
    del fields
    selected = feature_query(RNA, DOT_PLOT_FEATURES)
    return rna_view([
        placed(text_box("About this page", Q2_TEXT), [0, 0], [1100, 200]),
        placed(
            cell_filter("Choose groups to compare", [DIAGNOSIS, SITE, CELL_TYPE]),
            [1105, 0],
            [485, 200],
        ),
        placed(dot_plot("Selected genes by Diagnosis", DIAGNOSIS, selected), [0, 205], [500, 655]),
        placed(dot_plot("Selected genes by Biopsy Site", SITE, selected), [505, 205], [540, 325]),
        placed(
            dot_plot("Selected genes by Inflammation Status", INFLAMMATION, selected),
            [505, 535],
            [540, 325],
        ),
        placed(
            dot_plot(
                "Selected genes by Patient",
                PATIENT,
                selected,
                y_label_size=110,
                only=(PATIENT, Q2_PATIENTS),
            ),
            [1050, 205],
            [540, 655],
        ),
    ])


def build_q6_view(fields: dict[str, str]) -> dict[str, Any]:
    selected = feature_query(RNA, DOT_PLOT_FEATURES)
    per_disease = [
        placed(
            dot_plot(
                f"Selected genes in {diagnosis}",
                CELL_TYPE,
                selected,
                only=(DIAGNOSIS, [diagnosis]),
            ),
            [795 + index * 400, 505],
            [395, 355],
        )
        for index, diagnosis in enumerate(Q6_DIAGNOSES)
    ]
    return rna_view([
        placed(text_box("About this page", Q6_TEXT), [0, 0], [1100, 150]),
        placed(cell_filter("Focus on cell types", [CELL_TYPE]), [1105, 0], [485, 150]),
        placed(
            heat_map(
                "Marker genes by Diagnosis (average expression)",
                DIAGNOSIS,
                [fields[gene] for gene in MARKER_GENES],
            ),
            [0, 155],
            [790, 705],
        ),
        placed(
            box_plot(f"{Q6_BOX_GENE} by Diagnosis", DIAGNOSIS, fields[Q6_BOX_GENE]),
            [795, 155],
            [795, 345],
        ),
        *per_disease,
    ])


def build_q7_view(fields: dict[str, str]) -> dict[str, Any]:
    selected = feature_query(
        LinkedFeatures(RNA_DATASOURCE, SUBGROUP, COLLECTION_GENES),
        max(len(genes) for genes in GENE_COLLECTIONS.values()),
    )
    per_disease = [
        placed(
            dot_plot(
                f"Collection genes in {diagnosis}",
                CELL_TYPE,
                selected,
                only=(DIAGNOSIS, [diagnosis]),
            ),
            [795 + index * 400, 205],
            [395, 655],
        )
        for index, diagnosis in enumerate(Q6_DIAGNOSES)
    ]
    return rna_view(
        [
            placed(text_box("About this page", Q7_TEXT), [0, 0], [1100, 190]),
            placed(cell_filter("Focus on cell types", [CELL_TYPE]), [1105, 0], [485, 190]),
            placed(
                heat_map(
                    "Collection genes by Diagnosis (average expression)",
                    DIAGNOSIS,
                    [fields[gene] for gene in COLLECTION_GENES],
                    y_label_size=150,
                ),
                [0, 205],
                [790, 655],
            ),
            *per_disease,
        ],
        marker_panel(COLLECTION_GENES),
    )


def build_q8_view(fields: dict[str, str]) -> dict[str, Any]:
    selected = feature_query(RNA, DOT_PLOT_FEATURES)
    panels = [
        placed(
            dot_plot(
                Q8_GROUP_TITLES[group],
                CELL_GROUP,
                selected,
                y_label_size=110,
                only=(TNF_RESPONSE, [group]),
            ),
            [index * 795, 200],
            [790, 390],
        )
        for index, group in enumerate(Q8_GROUPS)
    ]
    return rna_view([
        placed(text_box("About this page", Q8_TEXT), [0, 0], [1100, 190]),
        placed(
            cell_filter(
                "Choose patients to compare",
                [DIAGNOSIS, HAD_TNF, TNF_RESPONSE],
                {DIAGNOSIS: ["RA"], HAD_TNF: ["Yes"], TNF_RESPONSE: Q8_GROUPS},
            ),
            [1105, 0],
            [485, 190],
        ),
        *panels,
        placed(
            composition_bar("Cell Group mix by TNF Treatment Response", TNF_RESPONSE, CELL_GROUP),
            [0, 600],
            [790, 260],
        ),
        placed(
            box_plot(
                f"{Q8_BOX_GENE} by TNF Treatment Response",
                TNF_RESPONSE,
                fields[Q8_BOX_GENE],
            ),
            [795, 600],
            [795, 260],
        ),
    ])


def build_q9_view(fields: dict[str, str]) -> dict[str, Any]:
    del fields
    selected = feature_query(RNA, DOT_PLOT_FEATURES)
    maps = [
        placed(
            umap(f"Cell Group: {Q8_GROUP_TITLES[group]}", CELL_GROUP, only=(TNF_RESPONSE, group)),
            [index * 530, 185],
            [520, 370],
        )
        for index, group in enumerate(Q8_GROUPS)
    ]
    plots = [
        placed(
            dot_plot(
                f"Selected genes: {Q8_GROUP_TITLES[group]}",
                CELL_GROUP,
                selected,
                y_label_size=110,
                only=(TNF_RESPONSE, [group]),
            ),
            [index * 795, 565],
            [790, 295],
        )
        for index, group in enumerate(Q8_GROUPS)
    ]
    return rna_view([
        placed(text_box("About this page", Q9_TEXT), [0, 0], [1100, 175]),
        placed(
            cell_filter(
                "Choose patients to compare",
                [DIAGNOSIS, HAD_TNF, TNF_RESPONSE],
                {DIAGNOSIS: ["RA"], HAD_TNF: ["Yes"], TNF_RESPONSE: Q8_GROUPS},
            ),
            [1105, 0],
            [485, 175],
        ),
        *maps,
        placed(
            composition_bar("Cell Group mix by TNF Treatment Response", TNF_RESPONSE, CELL_GROUP),
            [1060, 185],
            [530, 370],
        ),
        *plots,
    ])


VIEW_BUILDERS = {
    Q1_VIEW_NAME: build_q1_view,
    Q2_VIEW_NAME: build_q2_view,
    Q6_VIEW_NAME: build_q6_view,
    Q7_VIEW_NAME: build_q7_view,
    Q8_VIEW_NAME: build_q8_view,
    Q9_VIEW_NAME: build_q9_view,
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
        HAD_TNF, TNF_RESPONSE,
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
    absent = [gene for gene in PLOTTED_GENES.names if gene not in gene_names]
    if absent:
        raise ValueError(f"Genes missing from '{RNA_DATASOURCE}': {absent}")
    for field, categories in (
        (DIAGNOSIS, ["RA", "PsA"]),
        (HAD_TNF, ["Yes"]),
        (TNF_RESPONSE, Q8_GROUPS),
        (PATIENT, Q2_PATIENTS),
    ):
        values = next(column["values"] for column in cells["columns"] if column.get("field") == field)
        unknown = [category for category in categories if category not in values]
        if unknown:
            raise ValueError(f"'{field}' has no categories {unknown}")


def retitle_text(value: str) -> str:
    if value in EXACT_TITLES:
        return EXACT_TITLES[value]
    updated = value
    for field, title in sorted(TITLE_BY_FIELD.items(), key=lambda item: len(item[0]), reverse=True):
        updated = updated.replace(field, title)
    return updated


def retitle_existing_views(views: dict[str, Any]) -> int:
    """Point stored titles at the sheet's final names. Field ids in param stay put."""
    changed = 0
    for view in views.values():
        if not isinstance(view, dict):
            continue
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


def h5_stamp(project_path: Path) -> tuple[int, int]:
    stat = (project_path / "datafile.h5").stat()
    return stat.st_size, stat.st_mtime_ns


def add_synovial_question_views(project_path: Path, *, dry_run: bool = False) -> dict[str, Any]:
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
        fields = resolve_feature_fields(handle, datasources, PLOTTED_GENES)
    after_read = h5_stamp(project_path)
    if after_read != before:
        raise RuntimeError(f"Reading datafile.h5 changed it: {before} -> {after_read}")

    views = json.loads(views_path.read_text(encoding="utf-8"))
    state = json.loads(state_path.read_text(encoding="utf-8"))
    retitled = retitle_existing_views(views)
    all_views = list(state.get("all_views") or [])
    for name, build in VIEW_BUILDERS.items():
        views[name] = build(fields)
        if name not in all_views:
            all_views.append(name)
    state["all_views"] = all_views

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


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", type=Path, default=DEFAULT_PROJECT)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    print(json.dumps(add_synovial_question_views(args.project, dry_run=args.dry_run), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
