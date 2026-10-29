#!/usr/bin/env python3
"""Curate an MDV project from the Data Cleaning sheet and add question views.

The default is a dry run: it plans renames, soft-deletes, a views backup, and
question views, then writes a JSON report beside ``--output``. It does not
create the output project, edit the input project, or open datafile.h5 for
writing. ``--write`` copies the project to ``--output`` and applies the plan
there. datafile.h5 is never modified; size and mtime are checked before and
after.

Question views use columns and gene rows that already exist. Protein
concordance is added only when a protein table is already linked. Spatial
localisation, collection-score columns, and DEG columns are left out.
"""

from __future__ import annotations

import argparse
import copy
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import h5py

from mdvtools.examples.rename_softdelete_mdv_project_fields import (
    apply_plan,
    apply_soft_delete_plan,
    backup_project_jsons,
    build_plan,
    build_report,
    build_soft_delete_plan,
    load_views,
    print_report,
    resolve_datasources,
)
from mdvtools.template_views.create_question_views import (
    GENE_COLLECTIONS,
    GENES,
    NAME_FIELD,
    Q1_VIEW_NAME,
    Q2_VIEW_NAME,
    Q3_MATCHING,
    Q3_VIEW_NAME,
    Q6_VIEW_NAME,
    Q7_VIEW_NAME,
    LinkedFeatures,
    atomic_write_json,
    box_plot,
    cell_filter,
    composition_bar,
    compound_dot_plot,
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
from mdvtools.template_views.create_synovial_question_views import MARKER_GENES

Q8_VIEW_NAME = "Contrasting responders and non-responders"
Q9_VIEW_NAME = "Comparing groups: is the UMAP different between responders and non-responders?"

LAYOUT_WIDTH = 1590
DOT_PLOT_FEATURES = 5
PREFERRED_DIAGNOSES = ["Crohns disease", "Ulcerative colitis", "RA", "PsA"]
PREFERRED_RESPONSE_GROUPS = (
    ["responder", "non_responder"],
    ["R", "NR"],
    ["responsive", "non-responsive"],
)
ROLE_BY_NAME = {
    "cell group": "cell_group",
    "cell compartment": "cell_group",
    "cell type": "cell_type",
    "cell state": "cell_type",
    "cell subgroup": "cell_subgroup",
    "diagnosis": "diagnosis",
    "tissue": "tissue",
    "biopsy site": "tissue",
    "inflammation": "inflammation",
    "inflammation status": "inflammation",
    "tnf treatment response": "response",
    "responder (r/nr)": "response",
    "response": "response",
    "prior tnf inhibitor": "prior_tnf",
    "umap 1": "umap_x",
    "umap 2": "umap_y",
}
UMAP_FIELD_PAIRS = (
    ("rna:X_umap_mindist_0.25_1", "rna:X_umap_mindist_0.25_2"),
    ("X_umap_mindist_0.25_1", "X_umap_mindist_0.25_2"),
    ("X_umap_1", "X_umap_2"),
)


def h5_stamp(project_path: Path) -> tuple[int, int] | None:
    path = project_path / "datafile.h5"
    if not path.is_file():
        return None
    stat = path.stat()
    return stat.st_size, stat.st_mtime_ns


def require_same_h5(before: tuple[int, int] | None, after: tuple[int, int] | None, label: str) -> None:
    if before != after:
        raise RuntimeError(f"{label} datafile.h5 changed: {before} -> {after}")


def report_path_for(project: Path, output: Path, *, dry_run: bool) -> Path:
    suffix = "dry_run_report.json" if dry_run else "report.json"
    path = (output.parent / f"{output.name}_{suffix}").resolve()
    project_path = project.resolve()
    if path == project_path or project_path in path.parents:
        raise ValueError(f"Refusing to write the report inside the input project: {path}")
    if path == output.resolve() or output.resolve() in path.parents:
        raise ValueError(f"Refusing to write the report inside the output project: {path}")
    return path


def marker_names() -> list[str]:
    names: list[str] = []
    for name in [*MARKER_GENES, *GENES.names]:
        if name not in names:
            names.append(name)
    return names


def collection_names() -> list[str]:
    names: list[str] = []
    for genes in GENE_COLLECTIONS.values():
        for name in genes:
            if name not in names:
                names.append(name)
    return names


def text_column_values(column: dict[str, Any]) -> list[str]:
    values = column.get("values")
    if not isinstance(values, list):
        return []
    return [str(value) for value in values]


def pick_subgroup(subgroups: dict[str, Any]) -> str | None:
    keys = [key for key in subgroups if isinstance(key, str)]
    for hint in ("logged", "log", "expr"):
        for key in keys:
            if hint in key.casefold():
                return key
    return keys[0] if keys else None


def find_cells(datasources: list[dict[str, Any]]) -> dict[str, Any] | None:
    named = next((item for item in datasources if item.get("name") == "cells"), None)
    if named is not None:
        return named
    for item in datasources:
        for link in (item.get("links") or {}).values():
            subgroups = (link.get("rows_as_columns") or {}).get("subgroups") or {}
            if subgroups:
                return item
    return None


def linked_tables(
    cells: dict[str, Any], datasources: list[dict[str, Any]]
) -> list[tuple[str, str, dict[str, Any]]]:
    found: list[tuple[str, str, dict[str, Any]]] = []
    for name, link in (cells.get("links") or {}).items():
        subgroups = (link.get("rows_as_columns") or {}).get("subgroups") or {}
        subgroup = pick_subgroup(subgroups)
        table = next((item for item in datasources if item.get("name") == name), None)
        if subgroup and table is not None:
            fields = {column.get("field") for column in table.get("columns", [])}
            if NAME_FIELD in fields:
                found.append((name, subgroup, table))
    return found


def table_kind(name: str) -> int:
    folded = name.casefold()
    if "protein" in folded or folded.startswith("adt"):
        return 2
    if folded in {"rna", "genes", "gene"} or "gene" in folded:
        return 0
    return 1


def choose_umap(columns: list[dict[str, Any]]) -> tuple[str | None, str | None]:
    grouped: dict[str, list[dict[str, Any]]] = {"umap_x": [], "umap_y": []}
    fields = {str(column.get("field")): column for column in columns}
    for column in columns:
        role = column_role(column)
        if role in grouped:
            grouped[role].append(column)

    def prefer(candidates: list[dict[str, Any]], suffix: str) -> dict[str, Any] | None:
        if not candidates:
            return None
        for hint in ("mindist_0.25", "mindist"):
            for column in candidates:
                if hint in str(column.get("field")):
                    return column
        for column in candidates:
            if str(column.get("field")).endswith(suffix):
                return column
        return candidates[0]

    x_column = prefer(grouped["umap_x"], "_1")
    y_column = prefer(grouped["umap_y"], "_2")
    if x_column is not None and y_column is not None:
        return str(x_column["field"]), str(y_column["field"])
    for x_field, y_field in UMAP_FIELD_PAIRS:
        if x_field in fields and y_field in fields:
            return x_field, y_field
    return None, None


def column_role(column: dict[str, Any]) -> str | None:
    name = str(column.get("name") or "").strip().casefold()
    spaced = name.replace("_", " ")
    return ROLE_BY_NAME.get(name) or ROLE_BY_NAME.get(spaced)


def is_patient_column(column: dict[str, Any]) -> bool:
    name = str(column.get("name") or "").strip().casefold().replace("_", " ")
    field = str(column.get("field") or "").casefold()
    return name in {"patient", "patient id"} or field.endswith("patient_id") or field.endswith(":patient")


def choose_roles(columns: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    roles: dict[str, dict[str, Any]] = {}
    for column in columns:
        role = column_role(column)
        if role and role not in {"umap_x", "umap_y"} and role not in roles:
            roles[role] = column
    if "patient" not in roles:
        patient = next((column for column in columns if is_patient_column(column)), None)
        if patient is not None:
            roles["patient"] = patient
    return roles


def present_names(table: dict[str, Any], wanted: list[str]) -> list[str]:
    column = next((item for item in table.get("columns", []) if item.get("field") == NAME_FIELD), None)
    if column is None:
        return []
    available = set(text_column_values(column))
    return [name for name in wanted if name in available]


def preferred_values(values: list[str], preferred: list[str], limit: int) -> list[str]:
    chosen = [value for value in preferred if value in values]
    if chosen:
        return chosen[:limit]
    skipped = {"", "nd", "na", "nan", "none", "null"}
    return [value for value in values if value.strip().casefold() not in skipped][:limit]


def response_groups(values: list[str]) -> list[str]:
    folded = {value.casefold(): value for value in values}
    for groups in PREFERRED_RESPONSE_GROUPS:
        if all(group.casefold() in folded for group in groups):
            return [folded[group.casefold()] for group in groups]
    return preferred_values(values, [], 2)


def apply_planned_names(
    datasources: list[dict[str, Any]], outcomes: dict[str, list[dict[str, Any]]]
) -> list[dict[str, Any]]:
    updated = copy.deepcopy(datasources)
    by_name = {item.get("name"): item for item in updated}
    for datasource, results in outcomes.items():
        metadata = by_name.get(datasource)
        if metadata is None:
            continue
        columns = {column.get("field"): column for column in metadata.get("columns", [])}
        for result in results:
            column = columns.get(result.get("field"))
            if column is None:
                continue
            if result.get("action") == "rename" and result.get("status") in {"planned", "renamed"}:
                column["name"] = result["target_name"]
            if result.get("action") == "soft_delete" and result.get("status") in {"planned", "soft_deleted"}:
                column["deleted"] = True
    return updated


def surviving_columns(metadata: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        column
        for column in metadata.get("columns", [])
        if column.get("field") and not column.get("deleted")
    ]


def name_replacements(outcomes: dict[str, list[dict[str, Any]]]) -> list[tuple[str, str]]:
    grouped: dict[str, set[str]] = {}
    for results in outcomes.values():
        for item in results:
            if item.get("action") != "rename" or item.get("status") not in {"planned", "renamed"}:
                continue
            target = str(item.get("target_name") or "")
            for old in (item.get("current_name"), item.get("display_name"), item.get("field")):
                if isinstance(old, str) and old and old != target:
                    grouped.setdefault(old, set()).add(target)
    return [(old, next(iter(targets))) for old, targets in grouped.items() if len(targets) == 1]


def retitle_text(value: str, exact: dict[str, str], phrases: list[tuple[str, str]]) -> str:
    if value in exact:
        return exact[value]
    stripped = value.strip()
    if stripped in exact and stripped != value:
        return exact[stripped]
    updated = value
    for old, new in phrases:
        updated = updated.replace(old, new)
    return updated


def retitle_views(views: dict[str, Any], replacements: list[tuple[str, str]]) -> int:
    """Update stored titles and axis labels. Chart param values stay as field ids."""
    exact = dict(replacements)
    phrases = [
        pair for pair in replacements if len(pair[0]) >= 8 or ":" in pair[0]
    ]
    phrases.sort(key=lambda pair: len(pair[0]), reverse=True)
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
                    updated = retitle_text(title, exact, phrases)
                    if updated != title:
                        chart["title"] = updated
                        changed += 1
                axis = chart.get("axis")
                if not isinstance(axis, dict):
                    continue
                for spec in axis.values():
                    if isinstance(spec, dict) and isinstance(spec.get("label"), str):
                        updated = retitle_text(spec["label"], exact, phrases)
                        if updated != spec["label"]:
                            spec["label"] = updated
                            changed += 1
    return changed


def umap_chart(
    title: str,
    x_field: str,
    y_field: str,
    color_by: str,
    tooltip: str,
    only: tuple[str, str] | None = None,
) -> dict[str, Any]:
    chart = scatter(title, [x_field, y_field], color_by, tooltip)
    if only:
        chart["category_filters"] = [{"column": only[0], "category": [only[1]]}]
    return chart


def gene_panel(datasource: str, names: list[str]) -> list[dict[str, Any]]:
    del datasource
    return [
        placed(feature_list("Marker genes", names), [0, 0], [270, 400]),
        placed(feature_table("Click a gene to plot it"), [0, 405], [270, 455]),
    ]


def question_view(
    cells_name: str,
    gene_name: str,
    cells: list[dict[str, Any]],
    panel: list[dict[str, Any]],
) -> dict[str, Any]:
    return {
        "dataSources": {cells_name: {"panelWidth": 85}, gene_name: {"panelWidth": 15}},
        "initialCharts": {cells_name: cells, gene_name: panel},
    }


def column_label(column: dict[str, Any] | None, fallback: str) -> str:
    if column is None:
        return fallback
    name = str(column.get("name") or "").strip()
    return name or fallback


def plan_question_views(
    project: Path,
    datasources: list[dict[str, Any]],
    views: dict[str, Any],
    marker_genes: list[str] | None = None,
) -> dict[str, Any]:
    """Build question views from columns already in the project. The h5 is read-only."""
    skipped: list[dict[str, str]] = [
        {"question": "spatial localisation", "reason": "spatial charts are omitted"},
        {"question": "collection scores", "reason": "score columns are not created"},
        {"question": "DEG", "reason": "DEG columns are not created"},
    ]
    cells = find_cells(datasources)
    if cells is None:
        skipped.append({"question": "all", "reason": "no cell datasource with a gene link"})
        return {"added": [], "skipped": skipped, "views": {}, "retitled_charts": 0}

    tables = linked_tables(cells, datasources)
    gene_tables = [item for item in tables if table_kind(item[0]) < 2]
    protein_tables = [item for item in tables if table_kind(item[0]) == 2]
    if not gene_tables:
        skipped.append({"question": "all", "reason": "no rows-as-columns link to a gene table"})
        return {"added": [], "skipped": skipped, "views": {}, "retitled_charts": 0}
    gene_name, subgroup, gene_table = sorted(gene_tables, key=lambda item: table_kind(item[0]))[0]
    columns = surviving_columns(cells)
    roles = choose_roles(columns)
    umap_x, umap_y = choose_umap(columns)
    requested_markers = marker_names() if marker_genes is None else marker_genes
    markers = present_names(gene_table, requested_markers)
    collections = present_names(gene_table, collection_names())
    group = roles.get("cell_group") or roles.get("cell_subgroup") or roles.get("cell_type")
    dot_category = roles.get("cell_type") or roles.get("cell_subgroup") or roles.get("cell_group")
    built: dict[str, dict[str, Any]] = {}

    if not markers:
        skipped.append({"question": "gene questions", "reason": "none of the marker genes are in the gene table"})
        return {"added": [], "skipped": skipped, "views": {}, "retitled_charts": 0}

    before = h5_stamp(project)
    features = LinkedFeatures(gene_name, subgroup, list(dict.fromkeys([*markers, *collections])))
    try:
        with h5py.File(project / "datafile.h5", "r") as handle:
            fields = resolve_feature_fields(handle, datasources, features)
    except (OSError, ValueError, KeyError, IndexError) as exc:
        skipped.append({"question": "gene questions", "reason": f"could not read gene rows: {exc}"})
        return {"added": [], "skipped": skipped, "views": {}}
    require_same_h5(before, h5_stamp(project), "Reading")

    cells_name = str(cells.get("name"))
    group_field = str(group.get("field")) if group else ""
    group_label = column_label(group, "Cell Group")
    dot_field = str(dot_category.get("field")) if dot_category else ""
    dot_label = column_label(dot_category, group_label)
    selected = feature_query(LinkedFeatures(gene_name, subgroup, markers), DOT_PLOT_FEATURES)
    highlight = [gene for gene in ("CD14", "S100A8", "TNF") if gene in fields] or markers[:3]

    def accept(name: str, view: dict[str, Any]) -> None:
        if name in views:
            skipped.append({"question": name, "reason": "already in the project"})
            return
        built[name] = view

    if umap_x and umap_y and group_field:
        maps = [
            placed(
                umap_chart(gene, umap_x, umap_y, fields[gene], dot_field or group_field),
                [390 + index * 400, 125],
                [395, 340],
            )
            for index, gene in enumerate(highlight)
        ]
        accept(
            Q1_VIEW_NAME,
            question_view(
                cells_name,
                gene_name,
                [
                    placed(
                        text_box(
                            "About this page",
                            f"**Where is my gene of interest expressed across cell populations?** "
                            f"The first map shows **{group_label}**. The next maps show where "
                            f"{', '.join(highlight)} are expressed (darker = higher). In the **dot plot**, "
                            f"a bigger dot means more cells express the gene, and a darker dot means a "
                            f"higher average level. The **heat map** shows the average level of the marker "
                            f"genes in each {dot_label}.\n\n"
                            f"**To look at other genes:** click a gene in the list on the right "
                            f"(Ctrl/Cmd-click to pick several).",
                        ),
                        [0, 0],
                        [LAYOUT_WIDTH, 120],
                    ),
                    placed(umap_chart(group_label, umap_x, umap_y, group_field, dot_field or group_field), [0, 125], [385, 340]),
                    *maps,
                    placed(
                        dot_plot(f"Selected genes by {dot_label}", dot_field, selected),
                        [0, 475],
                        [700, 390],
                    ),
                    placed(
                        heat_map(
                            f"Marker genes by {group_label} (average expression)",
                            group_field,
                            [fields[gene] for gene in markers],
                        ),
                        [710, 475],
                        [880, 390],
                    ),
                ],
                gene_panel(gene_name, markers),
            ),
        )
    else:
        skipped.append({"question": Q1_VIEW_NAME, "reason": "needs a UMAP pair, a grouping column, and marker genes"})

    comparisons = [
        roles[role]
        for role in ("diagnosis", "tissue", "inflammation", "patient", "cell_type")
        if role in roles
    ]
    if comparisons and dot_field:
        comparison_fields = [str(column["field"]) for column in comparisons]
        charts = [
            placed(
                text_box(
                    "About this page",
                    "**Is the expression different between disease stages, patient groups, and tissue types?** "
                    "The **dot plots** compare the selected genes in the groups on this page. A bigger dot "
                    "means more cells express the gene, and a darker dot means a higher average level.\n\n"
                    "**To compare groups:** pick categories in **Choose groups to compare**. "
                    "**To change genes:** click genes in the list on the right.",
                ),
                [0, 0],
                [1100, 160],
            ),
            placed(cell_filter("Choose groups to compare", comparison_fields), [1105, 0], [485, 160]),
        ]
        for index, column in enumerate(comparisons[:4]):
            field = str(column["field"])
            label = column_label(column, field)
            only = None
            if roles.get("patient") is column:
                chosen = preferred_values(text_column_values(column), [], 5)
                if chosen:
                    only = (field, chosen)
            charts.append(
                placed(
                    dot_plot(f"Selected genes by {label}", field, selected, y_label_size=110, only=only),
                    [(index % 2) * 795, 170 + (index // 2) * 340],
                    [790, 330],
                )
            )
        accept(Q2_VIEW_NAME, question_view(cells_name, gene_name, charts, gene_panel(gene_name, markers)))
    else:
        skipped.append({"question": Q2_VIEW_NAME, "reason": "needs a grouping column and marker genes"})

    if protein_tables:
        protein_name, protein_subgroup, protein_table = protein_tables[0]
        pairs = [
            (gene, protein)
            for gene, protein in Q3_MATCHING
            if gene in fields and protein in set(present_names(protein_table, [protein]))
        ]
        protein_failed = False
        if pairs and dot_field:
            protein_features = LinkedFeatures(protein_name, protein_subgroup, [protein for _, protein in pairs])
            try:
                with h5py.File(project / "datafile.h5", "r") as handle:
                    resolve_feature_fields(handle, datasources, protein_features)
            except (OSError, ValueError, KeyError, IndexError) as exc:
                skipped.append({"question": Q3_VIEW_NAME, "reason": f"could not read protein rows: {exc}"})
                protein_failed = True
                pairs = []
            require_same_h5(before, h5_stamp(project), "Reading")
        if pairs and dot_field:
            accept(
                Q3_VIEW_NAME,
                {
                    "dataSources": {
                        cells_name: {"panelWidth": 70},
                        gene_name: {"panelWidth": 15},
                        protein_name: {"panelWidth": 15},
                    },
                    "initialCharts": {
                        cells_name: [
                            placed(
                                text_box(
                                    "About this page",
                                    "**Is the expression concordant with protein expression?** "
                                    "The left dot plot shows **RNA** and the right one shows **protein** "
                                    "for the pairs already measured in this project.",
                                ),
                                [0, 0],
                                [LAYOUT_WIDTH, 120],
                            ),
                            placed(
                                dot_plot(
                                    f"RNA by {dot_label}",
                                    dot_field,
                                    feature_query(
                                        LinkedFeatures(gene_name, subgroup, [gene for gene, _ in pairs]),
                                        DOT_PLOT_FEATURES,
                                    ),
                                ),
                                [0, 130],
                                [790, 500],
                            ),
                            placed(
                                dot_plot(
                                    f"Protein by {dot_label}",
                                    dot_field,
                                    feature_query(protein_features, DOT_PLOT_FEATURES),
                                ),
                                [800, 130],
                                [790, 500],
                            ),
                        ],
                        gene_name: gene_panel(gene_name, [gene for gene, _ in pairs]),
                        protein_name: [
                            placed(feature_list("Marker proteins", [protein for _, protein in pairs]), [0, 0], [270, 400]),
                            placed(feature_table("Click a protein to plot it"), [0, 405], [270, 455]),
                        ],
                    },
                },
            )
        elif not protein_failed:
            skipped.append({"question": Q3_VIEW_NAME, "reason": "protein table has none of the expected gene/protein pairs"})
    else:
        skipped.append({"question": Q3_VIEW_NAME, "reason": "no protein table"})

    diagnosis = roles.get("diagnosis")
    if diagnosis and dot_field:
        diagnosis_field = str(diagnosis["field"])
        diagnosis_values = preferred_values(text_column_values(diagnosis), PREFERRED_DIAGNOSES, 2)
        diagnosis_label = column_label(diagnosis, "Diagnosis")
        per_disease = placed(
            compound_dot_plot(
                f"Selected genes in {' and '.join(diagnosis_values)}",
                dot_field,
                selected,
                split=(diagnosis_field, diagnosis_values),
            ),
            [795, 505],
            [795, 355],
        )
        box_gene = "CD14" if "CD14" in fields else markers[0]
        accept(
            Q6_VIEW_NAME,
            question_view(
                cells_name,
                gene_name,
                [
                    placed(
                        text_box(
                            "About this page",
                            f"**How comparable is the expression across different diseases?** "
                            f"The **heat map** shows average marker expression by **{diagnosis_label}**. "
                            f"The **box plot** shows how {box_gene} varies within each group.",
                        ),
                        [0, 0],
                        [1100, 140],
                    ),
                    placed(cell_filter(f"Focus on {dot_label}", [dot_field]), [1105, 0], [485, 140]),
                    placed(
                        heat_map(
                            f"Marker genes by {diagnosis_label} (average expression)",
                            diagnosis_field,
                            [fields[gene] for gene in markers],
                        ),
                        [0, 150],
                        [790, 710],
                    ),
                    placed(
                        box_plot(f"{box_gene} by {diagnosis_label}", diagnosis_field, fields[box_gene]),
                        [795, 150],
                        [795, 345],
                    ),
                    per_disease,
                ],
                gene_panel(gene_name, markers),
            ),
        )
    else:
        skipped.append({"question": Q6_VIEW_NAME, "reason": "needs a diagnosis column and marker genes"})

    if diagnosis and len(collections) >= 2 and dot_field:
        diagnosis_field = str(diagnosis["field"])
        diagnosis_values = preferred_values(text_column_values(diagnosis), PREFERRED_DIAGNOSES, 2)
        diagnosis_label = column_label(diagnosis, "Diagnosis")
        collection_query = feature_query(
            LinkedFeatures(gene_name, subgroup, collections),
            len(collections),
        )
        per_disease = placed(
            compound_dot_plot(
                f"Collection genes in {' and '.join(diagnosis_values)}",
                dot_field,
                collection_query,
                split=(diagnosis_field, diagnosis_values),
            ),
            [800, 150],
            [790, 710],
        )
        accept(
            Q7_VIEW_NAME,
            question_view(
                cells_name,
                gene_name,
                [
                    placed(
                        text_box(
                            "About this page",
                            "**How comparable is the expression of gene collections across diseases?** "
                            "This page plots the **member genes** that are already in the project. "
                            "It does not add collection-score columns. A bigger dot means more cells "
                            "express the gene, and a darker dot means a higher average level.",
                        ),
                        [0, 0],
                        [1100, 140],
                    ),
                    placed(cell_filter(f"Focus on {dot_label}", [dot_field]), [1105, 0], [485, 140]),
                    placed(
                        heat_map(
                            f"Collection genes by {diagnosis_label} (average expression)",
                            diagnosis_field,
                            [fields[gene] for gene in collections],
                            y_label_size=150,
                        ),
                        [0, 150],
                        [790, 710],
                    ),
                    per_disease,
                ],
                gene_panel(gene_name, collections),
            ),
        )
    else:
        skipped.append({"question": Q7_VIEW_NAME, "reason": "needs a diagnosis column and at least two collection genes"})

    response = roles.get("response")
    groups = response_groups(text_column_values(response)) if response else []
    if response and groups and dot_field:
        response_field = str(response["field"])
        response_label = column_label(response, "Response")
        filters = [response_field]
        selected_filters: dict[str, list[str]] = {response_field: groups}
        if diagnosis:
            filters.insert(0, str(diagnosis["field"]))
            diagnosis_choice = preferred_values(text_column_values(diagnosis), PREFERRED_DIAGNOSES, 1)
            if diagnosis_choice:
                selected_filters[str(diagnosis["field"])] = diagnosis_choice
        prior = roles.get("prior_tnf")
        if prior and "Yes" in text_column_values(prior):
            filters.append(str(prior["field"]))
            selected_filters[str(prior["field"])] = ["Yes"]
        by_response = placed(
            compound_dot_plot(
                f"Selected genes by {response_label}",
                group_field or dot_field,
                selected,
                split=(response_field, [str(value) for value in groups[:2]]),
                y_label_size=110,
            ),
            [0, 200],
            [1590, 390],
        )
        box_gene = "TNF" if "TNF" in fields else markers[0]
        accept(
            Q8_VIEW_NAME,
            question_view(
                cells_name,
                gene_name,
                [
                    placed(
                        text_box(
                            "About this page",
                            f"**Contrasting response groups.** The **dot plot** shows the selected genes, "
                            f"with a panel for each value of **{response_label}**. The **bar** shows the mix of "
                            f"**{group_label}**, and the **box plot** shows how {box_gene} varies.",
                        ),
                        [0, 0],
                        [1100, 180],
                    ),
                    placed(
                        cell_filter("Choose patients to compare", filters, selected_filters),
                        [1105, 0],
                        [485, 180],
                    ),
                    by_response,
                    placed(
                        composition_bar(
                            f"{group_label} mix by {response_label}",
                            response_field,
                            group_field or dot_field,
                        ),
                        [0, 600],
                        [790, 260],
                    ),
                    placed(
                        box_plot(f"{box_gene} by {response_label}", response_field, fields[box_gene]),
                        [795, 600],
                        [795, 260],
                    ),
                ],
                gene_panel(gene_name, markers),
            ),
        )
    else:
        skipped.append({"question": Q8_VIEW_NAME, "reason": "needs an existing response column and marker genes"})

    if response and groups and umap_x and umap_y and group_field:
        response_field = str(response["field"])
        response_label = column_label(response, "Response")
        maps = [
            placed(
                umap_chart(
                    f"{group_label}: {group_value}",
                    umap_x,
                    umap_y,
                    group_field,
                    dot_field or group_field,
                    only=(response_field, group_value),
                ),
                [index * 530, 185],
                [520, 370],
            )
            for index, group_value in enumerate(groups[:2])
        ]
        accept(
            Q9_VIEW_NAME,
            question_view(
                cells_name,
                gene_name,
                [
                    placed(
                        text_box(
                            "About this page",
                            f"**Is the UMAP different between response groups?** The maps show each "
                            f"**{response_label}** group, coloured by **{group_label}**. A region that "
                            f"fills or empties points to a shift in cell composition. No DEG columns are added.",
                        ),
                        [0, 0],
                        [1100, 160],
                    ),
                    *maps,
                    placed(
                        composition_bar(
                            f"{group_label} mix by {response_label}",
                            response_field,
                            group_field,
                        ),
                        [1060, 185],
                        [530, 370],
                    ),
                ],
                gene_panel(gene_name, markers),
            ),
        )
    else:
        skipped.append({
            "question": Q9_VIEW_NAME,
            "reason": "needs a UMAP pair, a grouping column, and an existing response column",
        })

    return {
        "added": list(built),
        "skipped": skipped,
        "views": built,
        "gene_datasource": gene_name,
        "gene_subgroup": subgroup,
        "markers": markers,
        "collection_genes": collections,
    }


def backup_views(project: Path) -> Path:
    views_path = project / "views.json"
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    backup_path = views_path.with_name(f"{views_path.name}.bak.{timestamp}")
    backup_path.write_bytes(views_path.read_bytes())
    return backup_path


def planned_backup_paths(output: Path) -> dict[str, str]:
    return {
        "views_backup": str(output / "views.json.bak.<timestamp>"),
        "json_backup": str(output.parent / f"{output.name}_json_backup_<timestamp>"),
    }


def curate(
    project: Path,
    project_id: str,
    excel_path: Path,
    output: Path,
    *,
    write: bool,
) -> dict[str, Any]:
    project = project.expanduser().resolve()
    output = output.expanduser().resolve()
    excel_path = excel_path.expanduser().resolve()
    if not excel_path.is_file():
        raise FileNotFoundError(f"Spreadsheet does not exist: {excel_path}")
    if not project.is_dir():
        raise FileNotFoundError(f"Input project does not exist: {project}")
    if output == project or project in output.parents or output in project.parents:
        raise ValueError("The output project must be a different directory from the input project.")
    for name in ("datasources.json", "views.json", "state.json"):
        if not (project / name).is_file():
            raise FileNotFoundError(f"Not a complete MDV project; missing {project / name}")

    dry_run = not write
    source_before = h5_stamp(project)
    datasources = resolve_datasources(project, excel_path, project_id)
    views = load_views(project)
    warnings: dict[str, list[str]] = {}
    errors: dict[str, list[str]] = {}
    plans: dict[str, list[dict[str, Any]]] = {}
    kept_fields: dict[str, set[str]] = {}
    ambiguous_fields: dict[str, set[str]] = {}
    rename_outcomes: dict[str, list[dict[str, Any]]] = {}
    for datasource in datasources:
        plan, datasource_warnings, datasource_errors, kept, ambiguous = build_plan(
            project,
            datasource,
            excel_path,
            project_id,
            source_project_dir=project,
        )
        plans[datasource] = plan
        kept_fields[datasource] = kept
        ambiguous_fields[datasource] = ambiguous
        warnings[datasource] = datasource_warnings
        errors[datasource] = datasource_errors
        _summary, results, _actions = apply_plan(project, datasource, plan, True)
        rename_outcomes[datasource] = results

    source_datasources = json.loads((project / "datasources.json").read_text(encoding="utf-8"))
    planned_datasources = apply_planned_names(source_datasources, rename_outcomes)
    replacements = name_replacements(rename_outcomes)
    retitle_input = copy.deepcopy(views)
    retitled = retitle_views(retitle_input, replacements)
    question_plan = plan_question_views(project, planned_datasources, retitle_input)
    views_for_delete = copy.deepcopy(retitle_input)
    views_for_delete.update(question_plan["views"])

    outcomes: dict[str, list[dict[str, Any]]] = {}
    soft_delete_plans: dict[str, list[dict[str, Any]]] = {}
    for datasource in datasources:
        soft_delete_plans[datasource] = build_soft_delete_plan(
            project,
            datasource,
            kept_fields[datasource],
            ambiguous_fields[datasource],
            views_for_delete,
        )
        summary, _results, _actions = apply_plan(project, datasource, plans[datasource], True)
        delete_summary, delete_results, _delete_actions = apply_soft_delete_plan(
            project, datasource, soft_delete_plans[datasource], True
        )
        summary.update(delete_summary)
        outcomes[datasource] = rename_outcomes[datasource] + delete_results
        print(
            f"\nDatasource: {datasource}; dry_run=True; "
            f"renames={len(plans[datasource])}; "
            f"soft_deletes={len(soft_delete_plans[datasource])}; summary={dict(summary)}"
        )
    require_same_h5(source_before, h5_stamp(project), "Input")

    report = build_report(project, project_id, outcomes, warnings, errors)
    report["source_path"] = str(project)
    report["output_path"] = str(output)
    report["dry_run"] = dry_run
    report["h5_bytes"] = source_before[0] if source_before else None
    report["h5_mtime_ns"] = source_before[1] if source_before else None
    report["question_views"] = question_plan["added"]
    report["skipped_questions"] = question_plan["skipped"]
    report["retitled_charts"] = retitled
    report.update(planned_backup_paths(output))
    print_report(report)
    print("\nQuestion views")
    for name in question_plan["added"]:
        print(f"  add: {name}")
    for item in question_plan["skipped"]:
        print(f"  skip: {item['question']}: {item['reason']}")
    print(f"  retitle stored titles/labels: {retitled}")

    if dry_run:
        path = report_path_for(project, output, dry_run=True)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(report, indent=2), encoding="utf-8")
        print(f"\nSaved field report to: {path}")
        print("Dry run complete. The input project and datafile.h5 were not changed.")
        return report

    if output.exists():
        raise FileExistsError(f"Output project already exists: {output}")
    shutil.copytree(project, output)
    output_before = h5_stamp(output)
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    json_backup = output.parent / f"{output.name}_json_backup_{timestamp}"
    backup_project_jsons(output, json_backup)
    views_backup = backup_views(output)
    for datasource in datasources:
        apply_plan(output, datasource, plans[datasource], False)
        apply_soft_delete_plan(output, datasource, soft_delete_plans[datasource], False)

    output_views = load_views(output)
    retitle_views(output_views, replacements)
    for name, view in question_plan["views"].items():
        output_views[name] = view
    atomic_write_json(output / "views.json", output_views)
    state_path = output / "state.json"
    state = json.loads(state_path.read_text(encoding="utf-8"))
    all_views = list(state.get("all_views") or [])
    for name in output_views:
        if name not in all_views:
            all_views.append(name)
    state["all_views"] = all_views
    atomic_write_json(state_path, state)

    require_same_h5(output_before, h5_stamp(output), "Output")
    require_same_h5(source_before, h5_stamp(project), "Input")
    report["dry_run"] = False
    report["views_backup"] = str(views_backup)
    report["json_backup"] = str(json_backup)
    path = report_path_for(project, output, dry_run=False)
    path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"\nSaved field report to: {path}")
    print(f"Views backup: {views_backup}")
    print("Write complete. datafile.h5 was not modified.")
    return report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", type=Path, required=True, help="Input MDV project directory.")
    parser.add_argument("--project-id", required=True, help="Project id column in the Data Cleaning sheet.")
    parser.add_argument("--xlsx", type=Path, required=True, help="Workbook containing the Data Cleaning sheet.")
    parser.add_argument("--output", type=Path, required=True, help="Directory for the curated project.")
    parser.add_argument("--write", action="store_true", help="Copy the project to --output and apply the plan.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    curate(args.project, args.project_id, args.xlsx, args.output, write=args.write)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
