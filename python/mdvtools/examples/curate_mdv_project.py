#!/usr/bin/env python3
"""Curate an MDV project from the Data Cleaning sheet and add question views.

The default is a dry run: it plans renames, soft-deletes, a views backup, and
question views, then writes a JSON report beside ``--output``. It does not
create the output project, edit the input project, or open datafile.h5 for
writing. ``--write`` copies the project to ``--output`` and applies the plan
there. ``--no-copy`` applies the plan to ``--project`` instead. ``--link-heavy``
hard-links each top-level ``*.h5`` and symlinks ``spatial/`` and ``images/``
instead of copying those bytes. datafile.h5 is never modified; size and mtime
are checked before and after.

Question views use columns and gene rows that already exist. A protein
view is added only when a protein table is already linked, and it is placed
after the gene view. An existing gene-protein concordance view is resized and
limited to one matching pair. The CITE-seq concordance view and the
signature-collection view are left out. Spatial localisation,
collection-score columns, and DEG columns are left out.
"""

from __future__ import annotations

import argparse
import copy
import json
import os
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
    planned_soft_deletes,
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
from mdvtools.template_views.create_synovial_question_views import MARKER_GENES

Q8_VIEW_NAME = "Contrasting responders and non-responders"
Q8_TREATMENT_VIEW_NAME = "Contrasting responders and non-responders at baseline and post treatment"
Q9_VIEW_NAME = "Comparing groups: is the UMAP different between responders and non-responders?"
Q9_TREATMENT_VIEW_NAME = "Comparing groups: is the UMAP different between treatment states? DEG?"
PROTEIN_VIEW_NAME = "Where is my protein of interest expressed across cell populations?"
CONCORD_VIEW_NAME = "Is the gene expression concordant with protein expression?"
CELL_ATLAS_NAME = "Cell atlas"
QUESTION_ORDER = (
    Q1_VIEW_NAME,
    PROTEIN_VIEW_NAME,
    CONCORD_VIEW_NAME,
    Q2_VIEW_NAME,
    Q6_VIEW_NAME,
    Q8_TREATMENT_VIEW_NAME,
    Q8_VIEW_NAME,
    Q9_TREATMENT_VIEW_NAME,
    Q9_VIEW_NAME,
)
# These stay out of the project: the signature view has no collection column,
# and the CITE-seq view is a duplicate of CONCORD_VIEW_NAME.
OMITTED_VIEW_NAMES = (Q7_VIEW_NAME, Q3_VIEW_NAME)

# Chart sizes follow the hand-tuned Derm views on a laptop window about 1440px wide.
# One side list leaves the cells pane at 78%. Gene and protein lists together leave 66%.
DESIGN_WIDTH = 1440
CHART_GAP = 12
CELLS_SHARE = 78
CONCORD_SHARE = 66
SIDE_SHARE = 22
CONCORD_SIDE_SHARE = 17
TEXT_HEIGHT = 150
TEXT_TALL = 215
UMAP_HEIGHT = 400
DOT_HEIGHT = 620
DOT_Y_LABEL = 220
DOT_X_LABEL = 110
HEAT_HEIGHT = 430
FILTER_WIDTH = 420
FILTER_HEIGHT = 360
BAR_HEIGHT = 400
PROMPT_HEIGHT = 140
LIST_HEIGHT = 560
DETAIL_HEIGHT = 700
SPLIT_MAX = 5
OPENING_COMPARISON_GENES = 5
CONCORD_MAX_ITEMS = 1
COMPARE_ROLES = ("diagnosis", "condition", "treatment", "tissue", "inflammation")
PREFERRED_DIAGNOSES = ["Crohns disease", "Ulcerative colitis", "RA", "PsA"]
PREFERRED_RESPONSE_GROUPS = (
    ["responder", "non_responder"],
    ["R", "NR"],
    ["responsive", "non-responsive"],
    ["Remission", "Non_Remission"],
)
TREATMENT_STAGE_PAIRS = (
    ["Pre", "Post"],
    ["baseline", "post-treatment"],
)
ROLE_BY_NAME = {
    "cell group": "cell_group",
    "cell compartment": "cell_group",
    "bucket": "cell_group",
    "cell type": "cell_type",
    "cell state": "cell_type",
    "annotation": "cell_type",
    "cell subgroup": "cell_subgroup",
    "sub_bucket": "cell_subgroup",
    "diagnosis": "diagnosis",
    "diagnosed disease": "diagnosis",
    "disease": "diagnosis",
    "condition": "condition",
    "tissue": "tissue",
    "biopsy site": "tissue",
    "inflammation": "inflammation",
    "inflammation status": "inflammation",
    "presence of inflammation": "inflammation",
    "treatment": "treatment",
    "major cell type": "cell_type",
    "tnf treatment response": "response",
    "responder (r/nr)": "response",
    "response": "response",
    "remission status": "response",
    "prior tnf inhibitor": "prior_tnf",
    "umap 1": "umap_x",
    "umap 2": "umap_y",
    "umap dimension 1": "umap_x",
    "umap dimension 2": "umap_y",
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
    if folded in {"prot", "protein", "adt"} or "protein" in folded or folded.startswith("adt"):
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
    found: dict[str, list[dict[str, Any]]] = {}
    for column in columns:
        role = column_role(column)
        if role is None or role in {"umap_x", "umap_y"}:
            continue
        found.setdefault(role, []).append(column)
    roles: dict[str, dict[str, Any]] = {}
    for role, candidates in found.items():
        if role == "cell_type" and len(candidates) > 1:
            # Prefer a readable cell-type column over a fine cell-state list.
            roles[role] = min(candidates, key=_cell_type_rank)
        else:
            roles[role] = candidates[0]
    if "patient" not in roles:
        patient = next((column for column in columns if is_patient_column(column)), None)
        if patient is not None:
            roles["patient"] = patient
    return roles


def _cell_type_rank(column: dict[str, Any]) -> tuple[int, int]:
    count = len(text_column_values(column))
    return (count > SPLIT_MAX or count == 0, count or 999)


def opening_names(names: list[str], count: int) -> list[str]:
    """The names a list starts on. Charts follow whatever the user picks after that."""
    return names[:count]


def split_category_values(column: dict[str, Any]) -> list[str] | None:
    """Values for a split dot plot. One value is not split, and neither is a long list."""
    values = [value for value in text_column_values(column) if value.strip()]
    if 2 <= len(values) <= SPLIT_MAX:
        return values
    return None


def treatment_stages(column: dict[str, Any] | None) -> list[str]:
    """Baseline and post-treatment labels, when the column has that pair."""
    if column is None:
        return []
    folded = {value.casefold(): value for value in text_column_values(column)}
    for pair in TREATMENT_STAGE_PAIRS:
        if all(label.casefold() in folded for label in pair):
            return [folded[label.casefold()] for label in pair]
    return []


def content_width(share: int) -> int:
    return int(DESIGN_WIDTH * share / 100) - 8


def pair_width(share: int) -> int:
    return (content_width(share) - CHART_GAP) // 2


def present_names(table: dict[str, Any], wanted: list[str]) -> list[str]:
    column = next((item for item in table.get("columns", []) if item.get("field") == NAME_FIELD), None)
    if column is None:
        return []
    available = set(text_column_values(column))
    return [name for name in wanted if name in available]


def column_names(table: dict[str, Any]) -> set[str]:
    column = next((item for item in table.get("columns", []) if item.get("field") == NAME_FIELD), None)
    if column is None:
        return set()
    return set(text_column_values(column))


def stored_name(available: set[str], wanted: str) -> str | None:
    """Return the stored spelling, treating a hyphen and an underscore as the same."""
    folded = {name.casefold(): name for name in available}
    spellings = [wanted, wanted.replace("_", "-"), wanted.replace("-", "_")]
    for spelling in dict.fromkeys(spellings):
        found = folded.get(spelling.casefold())
        if found is not None:
            return found
    return None


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
    color_by: str | dict[str, Any],
    tooltip: str,
    only: tuple[str, str] | None = None,
) -> dict[str, Any]:
    """`color_by` is a field, or a `feature_query` for the selected gene (active link)."""
    chart = scatter(title, [x_field, y_field], color_by, tooltip)
    if only:
        chart["category_filters"] = [{"column": only[0], "category": [only[1]]}]
    return chart


def gene_panel(datasource: str, names: list[str], prompt: str | None = None) -> list[dict[str, Any]]:
    """`names` are the genes the list starts on; charts follow the genes picked in it.

    Question views leave `prompt` unset, so the list starts at the top of the pane.
    """
    del datasource
    return side_panel(names, prompt=prompt)


def side_panel(
    names: list[str],
    *,
    prompt: str | None = None,
    list_title: str = "Marker genes",
    kind: str = "gene",
    share: int = SIDE_SHARE,
) -> list[dict[str, Any]]:
    width = int(DESIGN_WIDTH * share / 100) - 16
    charts: list[dict[str, Any]] = []
    y = 0
    if prompt:
        charts.append(placed(text_box(prompt, prompt), [0, y], [width, PROMPT_HEIGHT]))
        y += PROMPT_HEIGHT + CHART_GAP
    charts.append(placed(feature_list(list_title, names), [0, y], [width, LIST_HEIGHT]))
    y += LIST_HEIGHT + CHART_GAP
    info_title = "Click a protein to plot it" if kind == "protein" else "Click a gene to plot it"
    charts.append(placed(gene_info(info_title), [0, y], [width, DETAIL_HEIGHT]))
    return charts


def wide_dot(title: str, category: str, features: dict[str, Any], **kwargs: Any) -> dict[str, Any]:
    """A full-width dot plot with room for category and gene names."""
    chart = compound_dot_plot(title, category, features, y_label_size=DOT_Y_LABEL, **kwargs)
    axis_x = chart["axis"]["x"]
    axis_x["size"] = DOT_X_LABEL
    axis_x["truncate_labels"] = False
    return chart


def pin_legends(view: dict[str, Any]) -> dict[str, Any]:
    """Colour legends sit bottom-right. Dot plots keep the default top legend."""
    for charts in view.get("initialCharts", {}).values():
        for chart in charts:
            legend = chart.get("color_legend")
            if chart.get("type") in {"dot_plot", "compound_dot_plot"}:
                if isinstance(legend, dict):
                    legend.pop("pos", None)
                continue
            size = chart.get("size")
            if isinstance(legend, dict) and legend.get("display") and isinstance(size, list) and len(size) == 2:
                width, height = size
                legend["pos"] = [max(8, int(width) - 180), max(8, int(height) - 120)]
    return view


def pin_project_legends(views: dict[str, Any]) -> None:
    """Move legends on every saved view, not only the newly built question views."""
    for view in views.values():
        if isinstance(view, dict):
            pin_legends(view)


def question_view(
    cells_name: str,
    gene_name: str,
    cells: list[dict[str, Any]],
    panel: list[dict[str, Any]],
    *,
    cells_share: int = CELLS_SHARE,
    side_share: int = SIDE_SHARE,
) -> dict[str, Any]:
    return {
        "dataSources": {
            cells_name: {"layout": "absolute", "panelWidth": cells_share},
            gene_name: {"layout": "absolute", "panelWidth": side_share},
        },
        "initialCharts": {cells_name: cells, gene_name: panel},
    }


class _Rows:
    """Place at most two charts on each row, inside the cells pane."""

    def __init__(self, share: int) -> None:
        self.share = share
        self.width = content_width(share)
        self.column = pair_width(share)
        self.y = 0
        self.charts: list[dict[str, Any]] = []

    def full(self, chart: dict[str, Any], height: int) -> None:
        self.charts.append(placed(chart, [0, self.y], [self.width, height]))
        self.y += height + CHART_GAP

    def pair(self, left: dict[str, Any], right: dict[str, Any] | None, height: int) -> None:
        self.charts.append(placed(left, [0, self.y], [self.column, height]))
        if right is not None:
            self.charts.append(placed(right, [self.column + CHART_GAP, self.y], [self.column, height]))
        self.y += height + CHART_GAP

    def selector(self, chart: dict[str, Any]) -> None:
        """Group selector after the other charts: narrower than the pane, and tall."""
        width = min(FILTER_WIDTH, self.width)
        self.charts.append(placed(chart, [0, self.y], [width, FILTER_HEIGHT]))
        self.y += FILTER_HEIGHT + CHART_GAP


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
        {"question": Q7_VIEW_NAME, "reason": "signature collection view is not created"},
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
    gene_features = LinkedFeatures(gene_name, subgroup, markers)
    one_gene = feature_query(gene_features, 1)
    many_genes = feature_query(gene_features, max(len(markers), 1))
    comparison_genes = opening_names(markers, OPENING_COMPARISON_GENES)
    counts = category_counts(datasources)
    y_field = group_field or dot_field
    axis_label = group_label if group_field else dot_label

    def accept(name: str, view: dict[str, Any]) -> None:
        built[name] = pin_legends(place_legends(view, counts))

    def split_chart(column: dict[str, Any], query: dict[str, Any]) -> dict[str, Any] | None:
        values = split_category_values(column)
        if values is None or not y_field:
            return None
        label = column_label(column, str(column.get("field")))
        return wide_dot(
            f"{axis_label}, split by {label}",
            y_field,
            query,
            split=(str(column["field"]), values),
        )

    comparisons = [roles[role] for role in COMPARE_ROLES if role in roles]

    if umap_x and umap_y and group_field and dot_field:
        rows = _Rows(CELLS_SHARE)
        rows.full(
            text_box(
                "About this page",
                f"**Where is my gene of interest expressed across cell populations?** "
                f"The first map shows **{group_label}**. The second map is the selected gene "
                f"(darker = higher). The **dot plot** and **heat map** show that one gene, by "
                f"**{dot_label}** and by **{group_label}**.\n\n"
                "**To look at another gene:** click one name in the gene list on the right.",
            ),
            TEXT_HEIGHT,
        )
        rows.pair(
            umap_chart(group_label, umap_x, umap_y, group_field, dot_field),
            umap_chart("Selected Gene Expression", umap_x, umap_y, one_gene, dot_field),
            UMAP_HEIGHT,
        )
        rows.full(heat_map(f"Selected gene by {group_label}", group_field, one_gene), HEAT_HEIGHT)
        rows.full(wide_dot(f"Selected gene by {dot_label}", dot_field, one_gene), DOT_HEIGHT)
        accept(
            Q1_VIEW_NAME,
            question_view(
                cells_name,
                gene_name,
                rows.charts,
                gene_panel(gene_name, [markers[0]]),
            ),
        )
    else:
        skipped.append({"question": Q1_VIEW_NAME, "reason": "needs a UMAP pair, a grouping column, and marker genes"})

    pairs: list[tuple[str, str]] = []
    protein_name = ""
    protein_subgroup = ""
    protein_failed = False
    if protein_tables and dot_field:
        protein_name, protein_subgroup, protein_table = protein_tables[0]
        pairs = []
        missing_genes: list[str] = []
        for gene, protein in Q3_MATCHING:
            stored = stored_name(column_names(protein_table), protein)
            if stored is None or not present_names(gene_table, [gene]):
                continue
            pairs.append((gene, stored))
            if gene not in fields:
                missing_genes.append(gene)
        if pairs:
            try:
                with h5py.File(project / "datafile.h5", "r") as handle:
                    if missing_genes:
                        fields.update(
                            resolve_feature_fields(
                                handle,
                                datasources,
                                LinkedFeatures(gene_name, subgroup, missing_genes),
                            )
                        )
                    resolve_feature_fields(
                        handle,
                        datasources,
                        LinkedFeatures(protein_name, protein_subgroup, [protein for _, protein in pairs]),
                    )
            except (OSError, ValueError, KeyError, IndexError) as exc:
                skipped.append({"question": CONCORD_VIEW_NAME, "reason": f"could not read protein rows: {exc}"})
                skipped.append({"question": PROTEIN_VIEW_NAME, "reason": f"could not read protein rows: {exc}"})
                protein_failed = True
                pairs = []
            else:
                pairs = [(gene, protein) for gene, protein in pairs if gene in fields]
            require_same_h5(before, h5_stamp(project), "Reading")
        if pairs and umap_x and umap_y and group_field:
            one_protein = feature_query(LinkedFeatures(protein_name, protein_subgroup, [pairs[0][1]]), 1)
            protein_rows = _Rows(CELLS_SHARE)
            protein_rows.full(
                text_box(
                    "About this page",
                    f"**Where is my protein of interest expressed across cell populations?** "
                    f"The first map shows **{group_label}**. The second map is the selected protein "
                    f"(darker = higher). The **dot plot** and **heat map** show that one protein, by "
                    f"**{dot_label}** and by **{group_label}**.\n\n"
                    "**To look at another protein:** click one name in the protein list on the right.",
                ),
                TEXT_HEIGHT,
            )
            protein_rows.pair(
                umap_chart(group_label, umap_x, umap_y, group_field, dot_field),
                umap_chart("Selected Protein Expression", umap_x, umap_y, one_protein, dot_field),
                UMAP_HEIGHT,
            )
            protein_rows.full(
                heat_map(f"Selected protein by {group_label}", group_field, one_protein),
                HEAT_HEIGHT,
            )
            protein_rows.full(
                wide_dot(f"Selected protein by {dot_label}", dot_field, one_protein),
                DOT_HEIGHT,
            )
            accept(
                PROTEIN_VIEW_NAME,
                question_view(
                    cells_name,
                    protein_name,
                    protein_rows.charts,
                    side_panel(
                        [pairs[0][1]],
                        list_title="Marker proteins",
                        kind="protein",
                    ),
                ),
            )
        elif pairs and not protein_failed:
            skipped.append({
                "question": PROTEIN_VIEW_NAME,
                "reason": "needs a UMAP pair and a grouping column",
            })
        elif not protein_failed:
            skipped.append({
                "question": PROTEIN_VIEW_NAME,
                "reason": "protein table has none of the expected gene/protein pairs",
            })
    elif protein_tables:
        skipped.append({"question": PROTEIN_VIEW_NAME, "reason": "needs a grouping column"})
    else:
        skipped.append({"question": PROTEIN_VIEW_NAME, "reason": "no protein table"})

    if comparisons and y_field:
        rows = _Rows(CELLS_SHARE)
        rows.full(
            text_box(
                "About this page",
                "**Is the expression different between disease stages, patient groups, and tissue types?** "
                "A category with several groups is a **split dot plot**: one panel per group, on one colour "
                "and size scale. A category with only one value stays a single panel. Cell type stays one "
                "panel when there are many types.\n\n"
                "The page starts on five genes. **To plot more:** click further genes in the list on the right.",
            ),
            TEXT_TALL,
        )
        filter_fields = [str(column["field"]) for column in comparisons]
        patient = roles.get("patient")
        if patient is not None:
            filter_fields.append(str(patient["field"]))
        response_column = roles.get("response")
        if response_column is not None:
            filter_fields.append(str(response_column["field"]))
        group_filter = cell_filter("Choose groups to compare", filter_fields)
        pending: list[dict[str, Any]] = []
        for column in comparisons:
            chart = split_chart(column, many_genes)
            if chart is None and len([value for value in text_column_values(column) if value.strip()]) == 1:
                label = column_label(column, str(column.get("field")))
                chart = wide_dot(
                    f"Selected genes by {label}",
                    str(column["field"]),
                    many_genes,
                )
            if chart is not None:
                pending.append(chart)
        if dot_field:
            pending.append(wide_dot(f"Selected genes by {dot_label}", dot_field, many_genes))
        for chart in pending:
            rows.full(chart, DOT_HEIGHT)
        rows.selector(group_filter)
        accept(Q2_VIEW_NAME, question_view(cells_name, gene_name, rows.charts, gene_panel(gene_name, comparison_genes)))
    else:
        skipped.append({"question": Q2_VIEW_NAME, "reason": "needs a grouping column and marker genes"})

    if pairs and dot_field and not protein_failed:
        concord = _Rows(CONCORD_SHARE)
        concord.full(
            text_box(
                "About this page",
                "**Is the gene expression concordant with protein expression?** "
                f"The charts start on one pair, **{pairs[0][0]}** and **{pairs[0][1]}**. "
                "The left side is RNA and the right side is protein. Similar patterns mean RNA and protein agree.\n\n"
                "**To compare another pair:** click that gene and its matching protein. "
                "The charts follow the one name selected in each list.",
            ),
            TEXT_HEIGHT,
        )
        focus_filter = cell_filter(f"Focus on {dot_label}", [dot_field])
        gene_query = feature_query(LinkedFeatures(gene_name, subgroup, [pairs[0][0]]), CONCORD_MAX_ITEMS)
        protein_query = feature_query(
            LinkedFeatures(protein_name, protein_subgroup, [pairs[0][1]]),
            CONCORD_MAX_ITEMS,
        )
        concord.full(wide_dot(f"RNA by {dot_label}", dot_field, gene_query), DOT_HEIGHT)
        concord.full(wide_dot(f"Protein by {dot_label}", dot_field, protein_query), DOT_HEIGHT)
        if umap_x and umap_y:
            concord.pair(
                umap_chart(
                    "Selected Gene Expression",
                    umap_x,
                    umap_y,
                    feature_query(LinkedFeatures(gene_name, subgroup, [pairs[0][0]]), 1),
                    dot_field,
                ),
                umap_chart(
                    "Selected Protein Expression",
                    umap_x,
                    umap_y,
                    feature_query(LinkedFeatures(protein_name, protein_subgroup, [pairs[0][1]]), 1),
                    dot_field,
                ),
                UMAP_HEIGHT,
            )
        concord.selector(focus_filter)
        accept(
            CONCORD_VIEW_NAME,
            {
                "dataSources": {
                    cells_name: {"layout": "absolute", "panelWidth": CONCORD_SHARE},
                    gene_name: {"layout": "absolute", "panelWidth": CONCORD_SIDE_SHARE},
                    protein_name: {"layout": "absolute", "panelWidth": CONCORD_SIDE_SHARE},
                },
                "initialCharts": {
                    cells_name: concord.charts,
                    gene_name: side_panel(
                        [pairs[0][0]],
                        share=CONCORD_SIDE_SHARE,
                    ),
                    protein_name: side_panel(
                        [pairs[0][1]],
                        list_title="Marker proteins",
                        kind="protein",
                        share=CONCORD_SIDE_SHARE,
                    ),
                },
            },
        )
    elif not protein_tables:
        skipped.append({"question": CONCORD_VIEW_NAME, "reason": "no protein table"})
    elif not protein_failed:
        skipped.append({"question": CONCORD_VIEW_NAME, "reason": "protein table has none of the expected gene/protein pairs"})

    split_columns = [column for column in comparisons if split_category_values(column)]
    if split_columns and y_field:
        rows = _Rows(CELLS_SHARE)
        rows.full(
            text_box(
                "About this page",
                "**How comparable is the expression across different diseases?** "
                "Each comparison is a **split dot plot** with a panel per group, and a **heat map** "
                "under it showing the average of the selected genes in that same grouping.\n\n"
                "The page starts on five genes. **To plot more:** click further genes in the list on the right.",
            ),
            TEXT_TALL,
        )
        group_filter = cell_filter(
            "Choose groups to compare", [str(column["field"]) for column in split_columns]
        )
        for column in split_columns:
            label = column_label(column, str(column.get("field")))
            chart = split_chart(column, many_genes)
            if chart is None:
                continue
            rows.full(chart, DOT_HEIGHT)
            rows.full(
                heat_map(f"Selected genes by {label}", str(column["field"]), many_genes),
                HEAT_HEIGHT,
            )
        rows.selector(group_filter)
        accept(Q6_VIEW_NAME, question_view(cells_name, gene_name, rows.charts, gene_panel(gene_name, comparison_genes)))
    else:
        skipped.append({"question": Q6_VIEW_NAME, "reason": "needs a diagnosis column and marker genes"})

    response = roles.get("response")
    treatment = roles.get("treatment")
    groups = response_groups(text_column_values(response)) if response else []
    stages = treatment_stages(treatment)
    if response and groups and y_field:
        response_field = str(response["field"])
        response_label = column_label(response, "Response")
        filters = [response_field]
        selected_filters: dict[str, list[str]] = {response_field: groups}
        if roles.get("diagnosis") is not None:
            filters.insert(0, str(roles["diagnosis"]["field"]))
        if treatment is not None:
            filters.append(str(treatment["field"]))
        rows = _Rows(CELLS_SHARE)
        if stages and treatment is not None:
            treatment_field = str(treatment["field"])
            treatment_label = column_label(treatment, "Treatment")
            rows.full(
                text_box(
                    "About this page",
                    f"**Contrasting {response_label} before and after treatment.** "
                    f"Each **dot plot** is limited to one **{treatment_label}** time and has a panel per "
                    f"**{response_label}** group, on a shared colour and size scale.\n\n"
                    "The page starts on five genes. **To plot more:** click further genes in the list on the right.",
                ),
                TEXT_TALL,
            )
            group_filter = cell_filter("Choose groups to compare", filters, selected_filters)
            stage_charts = []
            for stage in stages[:2]:
                chart = wide_dot(
                    f"{stage}: {response_label}",
                    y_field,
                    many_genes,
                    split=(response_field, [str(value) for value in groups[:2]]),
                )
                chart["category_filters"] = [{"column": treatment_field, "category": [stage]}]
                stage_charts.append(chart)
            for chart in stage_charts:
                rows.full(chart, DOT_HEIGHT)
            rows.selector(group_filter)
            q8_name = Q8_TREATMENT_VIEW_NAME
        else:
            rows.full(
                text_box(
                    "About this page",
                    f"**Contrasting response groups.** The **dot plot** shows the selected genes, with a panel "
                    f"for each value of **{response_label}**. The **bar** shows the mix of **{axis_label}**, "
                    f"and the **box plot** shows how {('TNF' if 'TNF' in fields else markers[0])} varies.\n\n"
                    "The page starts on five genes. **To plot more:** click further genes in the list on the right.",
                ),
                TEXT_TALL,
            )
            group_filter = cell_filter("Choose groups to compare", filters, selected_filters)
            rows.full(
                wide_dot(
                    f"{axis_label}, split by {response_label}",
                    y_field,
                    many_genes,
                    split=(response_field, [str(value) for value in groups[:2]]),
                ),
                DOT_HEIGHT,
            )
            box_gene = "TNF" if "TNF" in fields else markers[0]
            rows.pair(
                composition_bar(f"{axis_label} mix by {response_label}", response_field, y_field),
                box_plot(f"{box_gene} by {response_label}", response_field, fields[box_gene]),
                BAR_HEIGHT,
            )
            rows.selector(group_filter)
            q8_name = Q8_VIEW_NAME
        accept(q8_name, question_view(cells_name, gene_name, rows.charts, gene_panel(gene_name, comparison_genes)))
    else:
        skipped.append({"question": Q8_VIEW_NAME, "reason": "needs an existing response column and marker genes"})

    stage_column = treatment if stages else response if groups else None
    stage_values = stages if stages else groups[:2]
    if stage_column is not None and len(stage_values) >= 2 and umap_x and umap_y and group_field and y_field:
        stage_field = str(stage_column["field"])
        stage_label = column_label(stage_column, "Treatment" if stages else "Response")
        rows = _Rows(CELLS_SHARE)
        rows.full(
            text_box(
                "About this page",
                f"**Is the UMAP different between {stage_label} groups?** The two maps show "
                f"**{stage_values[0]}** and **{stage_values[1]}**, coloured by **{group_label}**. "
                f"The bar shows that mix. The **dot plot** has a panel for each group.\n\n"
                "**Differentially expressed genes** are not added, so there is no volcano plot. "
                "The page starts on five genes. **To plot more:** click further genes in the list on the right.",
            ),
            TEXT_TALL,
        )
        filter_fields = [stage_field]
        if response is not None and str(response["field"]) not in filter_fields:
            filter_fields.append(str(response["field"]))
        group_filter = cell_filter("Choose groups to compare", filter_fields)
        rows.pair(
            umap_chart(
                str(stage_values[0]),
                umap_x,
                umap_y,
                group_field,
                dot_field or group_field,
                only=(stage_field, str(stage_values[0])),
            ),
            umap_chart(
                str(stage_values[1]),
                umap_x,
                umap_y,
                group_field,
                dot_field or group_field,
                only=(stage_field, str(stage_values[1])),
            ),
            UMAP_HEIGHT,
        )
        box_gene = "TNF" if "TNF" in fields else markers[0]
        rows.pair(
            composition_bar(f"{group_label} mix by {stage_label}", stage_field, group_field),
            box_plot(f"{box_gene} by {stage_label}", stage_field, fields[box_gene]),
            BAR_HEIGHT,
        )
        rows.full(
            wide_dot(
                f"{axis_label}, split by {stage_label}",
                y_field,
                many_genes,
                split=(stage_field, [str(value) for value in stage_values[:2]]),
            ),
            DOT_HEIGHT,
        )
        rows.selector(group_filter)
        q9_name = Q9_TREATMENT_VIEW_NAME if stages else Q9_VIEW_NAME
        accept(q9_name, question_view(cells_name, gene_name, rows.charts, gene_panel(gene_name, comparison_genes)))
    else:
        skipped.append({
            "question": Q9_VIEW_NAME,
            "reason": "needs a UMAP pair, a grouping column, and an existing response or treatment column",
        })

    return {
        "added": [name for name in built if name not in views],
        "updated": [name for name in built if name in views],
        "removed": [name for name in OMITTED_VIEW_NAMES if name in views],
        "skipped": skipped,
        "views": built,
        "gene_datasource": gene_name,
        "gene_subgroup": subgroup,
        "markers": markers,
        "collection_genes": collections,
    }


def arrange_views(all_views: list[str], built: dict[str, Any]) -> list[str]:
    """Cell atlas, then the question views, then every other saved view."""
    names = [name for name in all_views if name not in OMITTED_VIEW_NAMES]
    questions = [name for name in QUESTION_ORDER if name in built]
    for name in built:
        if name not in questions:
            questions.append(name)
    others = [name for name in names if name not in questions and name != CELL_ATLAS_NAME]
    atlas = [CELL_ATLAS_NAME] if CELL_ATLAS_NAME in names else []
    return atlas + questions + others


def backup_views(project: Path) -> Path:
    views_path = project / "views.json"
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    backup_path = views_path.with_name(f"{views_path.name}.bak.{timestamp}")
    backup_path.write_bytes(views_path.read_bytes())
    return backup_path


HEAVY_DIR_NAMES = ("images", "spatial")


def ignore_unreadable(directory: str, names: list[str]) -> list[str]:
    """Names copytree should skip because this process cannot list the directory."""
    skipped: list[str] = []
    for name in names:
        path = Path(directory) / name
        if path.is_dir() and not os.access(path, os.R_OK | os.X_OK):
            print(f"Skipping unreadable directory: {path}")
            skipped.append(name)
    return skipped


def ignore_for_copy(root: Path, *, link_heavy: bool):
    """Skip unreadable directories, and with ``link_heavy`` the heavy root entries."""
    root = root.resolve()

    def ignore(directory: str, names: list[str]) -> list[str]:
        skipped = ignore_unreadable(directory, names)
        if not link_heavy or Path(directory).resolve() != root:
            return skipped
        for name in names:
            if name in skipped:
                continue
            path = Path(directory) / name
            if name in HEAVY_DIR_NAMES or path.suffix == ".h5":
                skipped.append(name)
        return skipped

    return ignore


def heavy_link_plan(project: Path) -> dict[str, list[str]]:
    """Top-level directories to symlink and h5 files to hard-link."""
    return {
        "symlink_dirs": [name for name in HEAVY_DIR_NAMES if (project / name).exists()],
        "hardlink_files": [path.name for path in sorted(project.glob("*.h5")) if path.is_file()],
    }


def link_heavy_paths(project: Path, output: Path) -> dict[str, list[str]]:
    """Share the heavy source files. Raises instead of falling back to a byte copy."""
    plan = heavy_link_plan(project)
    for name in plan["symlink_dirs"]:
        source = project / name
        destination = output / name
        destination.symlink_to(source)
        print(f"Symlinked {destination} -> {source}")
    for name in plan["hardlink_files"]:
        source = project / name
        destination = output / name
        try:
            os.link(source, destination)
        except OSError as exc:
            raise RuntimeError(
                f"Could not hard-link {source} to {destination}. The h5 file was not copied."
            ) from exc
        if destination.stat().st_ino != source.stat().st_ino:
            raise RuntimeError(f"Hard link {destination} does not share an inode with {source}.")
        print(f"Hard-linked {destination} -> {source}")
    return plan


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
    marker_genes: list[str] | None = None,
    link_heavy: bool = False,
    copy_project: bool = True,
) -> dict[str, Any]:
    project = project.expanduser().resolve()
    output = project if not copy_project else output.expanduser().resolve()
    excel_path = excel_path.expanduser().resolve()
    if not excel_path.is_file():
        raise FileNotFoundError(f"Spreadsheet does not exist: {excel_path}")
    if not project.is_dir():
        raise FileNotFoundError(f"Input project does not exist: {project}")
    if copy_project and (output == project or project in output.parents or output in project.parents):
        raise ValueError("The output project must be a different directory from the input project.")
    for name in ("datasources.json", "views.json", "state.json"):
        if not (project / name).is_file():
            raise FileNotFoundError(f"Not a complete MDV project; missing {project / name}")

    dry_run = not write
    source_before = h5_stamp(project)
    datasources, skipped_datasources = resolve_datasources(project, excel_path, project_id)
    if skipped_datasources:
        print(
            "Skipping datasources that are not in the project: "
            + ", ".join(skipped_datasources)
        )
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
        # Question views aren't planned yet, so this soft-delete list can be a
        # little too long; the final check below uses the real one.
        first_soft_deletes = build_soft_delete_plan(project, datasource, kept, ambiguous, views)
        _summary, results, _actions = apply_plan(
            project, datasource, plan, True, planned_soft_deletes(first_soft_deletes)
        )
        rename_outcomes[datasource] = results

    source_datasources = json.loads((project / "datasources.json").read_text(encoding="utf-8"))
    planned_datasources = apply_planned_names(source_datasources, rename_outcomes)
    replacements = name_replacements(rename_outcomes)
    retitle_input = copy.deepcopy(views)
    retitled = retitle_views(retitle_input, replacements)
    question_plan = plan_question_views(project, planned_datasources, retitle_input, marker_genes)
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
        summary, results, _actions = apply_plan(
            project,
            datasource,
            plans[datasource],
            True,
            planned_soft_deletes(soft_delete_plans[datasource]),
        )
        if [item["status"] for item in results] != [
            item["status"] for item in rename_outcomes[datasource]
        ]:
            print(
                f"Warning: {datasource}: a rename is blocked because a question view "
                "keeps a column with the same name visible."
            )
        delete_summary, delete_results, _delete_actions = apply_soft_delete_plan(
            project, datasource, soft_delete_plans[datasource], True
        )
        summary.update(delete_summary)
        outcomes[datasource] = results + delete_results
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
    report["skipped_datasources"] = skipped_datasources
    report["question_views"] = question_plan["added"]
    report["updated_question_views"] = question_plan.get("updated", [])
    report["skipped_questions"] = question_plan["skipped"]
    report["retitled_charts"] = retitled
    report.update(planned_backup_paths(output))
    if link_heavy:
        report["heavy_links"] = heavy_link_plan(project)
        print(
            "Heavy entries will be linked: "
            f"symlinks={report['heavy_links']['symlink_dirs']}; "
            f"hardlinks={report['heavy_links']['hardlink_files']}"
        )
    print_report(report)
    print("\nQuestion views")
    for name in question_plan["added"]:
        print(f"  add: {name}")
    for name in question_plan.get("updated", []):
        print(f"  update: {name}")
    for name in question_plan.get("removed", []):
        print(f"  remove: {name}")
    for item in question_plan["skipped"]:
        print(f"  skip: {item['question']}: {item['reason']}")
    print(f"  retitle stored titles/labels: {retitled}")
    if marker_genes is not None:
        found = set(question_plan.get("markers", []))
        missing = [name for name in marker_genes if name not in found]
        print(f"  markers found: {len(found)} of {len(marker_genes)}")
        if missing:
            print(f"  markers not in the gene table: {', '.join(missing)}")

    if dry_run:
        path = report_path_for(project, output, dry_run=True)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(report, indent=2), encoding="utf-8")
        print(f"\nSaved field report to: {path}")
        print("Dry run complete. The input project and datafile.h5 were not changed.")
        return report

    if copy_project:
        if output.exists():
            raise FileExistsError(f"Output project already exists: {output}")
        shutil.copytree(project, output, symlinks=True, ignore=ignore_for_copy(project, link_heavy=link_heavy))
        try:
            if link_heavy:
                link_heavy_paths(project, output)
        except Exception:
            shutil.rmtree(output)
            raise
    elif link_heavy:
        print("Keeping the project in place. Heavy files stay where they are.")
    output_before = h5_stamp(output)
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    json_backup = output.parent / f"{output.name}_json_backup_{timestamp}"
    backup_project_jsons(output, json_backup)
    views_backup = backup_views(output)
    for datasource in datasources:
        apply_plan(
            output,
            datasource,
            plans[datasource],
            False,
            planned_soft_deletes(soft_delete_plans[datasource]),
        )
        apply_soft_delete_plan(output, datasource, soft_delete_plans[datasource], False)

    output_views = load_views(output)
    retitle_views(output_views, replacements)
    for name in OMITTED_VIEW_NAMES:
        output_views.pop(name, None)
    for name, view in question_plan["views"].items():
        output_views[name] = view
    pin_project_legends(output_views)
    atomic_write_json(output / "views.json", output_views)
    state_path = output / "state.json"
    state = json.loads(state_path.read_text(encoding="utf-8"))
    all_views = [name for name in state.get("all_views") or [] if isinstance(name, str)]
    state["all_views"] = arrange_views(all_views, question_plan["views"])
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
    parser.add_argument(
        "--output",
        type=Path,
        help="Directory for the curated copy. Required unless --no-copy is set.",
    )
    parser.add_argument("--write", action="store_true", help="Apply the plan. Copies the project unless --no-copy is set.")
    parser.add_argument(
        "--no-copy",
        action="store_true",
        help="Apply the plan to --project instead of copying it to --output.",
    )
    parser.add_argument(
        "--link-heavy",
        action="store_true",
        help="With --write, hard-link top-level *.h5 files and symlink spatial/ and images/.",
    )
    parser.add_argument(
        "--markers",
        help='Comma-separated marker genes for the question views, e.g. "PTPRC, CD3E, CD4". '
        "Defaults to a general immune and stromal list.",
    )
    return parser.parse_args()


def parse_markers(text: str | None) -> list[str] | None:
    if text is None:
        return None
    names = [name.strip() for name in text.split(",") if name.strip()]
    return list(dict.fromkeys(names))


def main() -> int:
    args = parse_args()
    if args.no_copy:
        if args.output is not None and args.output.expanduser().resolve() != args.project.expanduser().resolve():
            raise SystemExit("--output must be omitted or equal to --project when --no-copy is set.")
        output = args.project
    elif args.output is None:
        raise SystemExit("--output is required unless --no-copy is set.")
    else:
        output = args.output
    curate(
        args.project,
        args.project_id,
        args.xlsx,
        output,
        write=args.write,
        marker_genes=parse_markers(args.markers),
        link_heavy=args.link_heavy,
        copy_project=not args.no_copy,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
