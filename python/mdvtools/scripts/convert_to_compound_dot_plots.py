#!/usr/bin/env python3
"""Merge groups of legacy dot plots into one compound dot plot, in views.json only.

A group is two or more `dot_plot` charts in the same view and datasource that plot the same
rows and fields (`param`), each restricted with a background filter to values of the same
column. Each group becomes one `compound_dot_plot` split by that column, with a panel per value,
in the area the group covered. A dot plot on its own (such as the Patient plot) is left alone.

Text boxes in a changed view get the wording the question-view scripts now use. Wording that
was edited by hand is left as it is and reported.

Stdlib only, so it runs with any Python 3.9+. Dry run unless --write; --write backs up
views.json first. datafile.h5, datasources.json and state.json are not changed.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

LEGACY = "dot_plot"
COMPOUND = "compound_dot_plot"
NO_COLUMN = "__none__"
# The compound plot's palettes; other legacy palettes fall back to the first.
PALETTES = ("reds", "blue_white_red")

# Settings carried over from the first plot in a group; the compound plot fills in the rest.
KEPT_KEYS = ("id", "legend", "param", "y_axis_order", "threshold_value")

# Old "About this page" wording -> new, for the question views made by
# create_question_views.py, create_synovial_question_views.py and curate_mdv_project.py.
TEXT_UPDATES: list[tuple[re.Pattern[str], str]] = [
    (
        re.compile(
            re.escape("**To change a dot plot's disease:** click its settings icon (⚙) and choose a "
                      "diagnosis under **Background Filter Categories**.")
        ),
        "**To compare other diseases:** click the dot plot's settings icon (⚙) and choose "
        "diagnoses under **Split values**.",
    ),
    (
        re.compile(r"The two \*\*dot plots\*\* show (.+?),? for ([^.:]+? and [^.:]+?)(?=[.:])"),
        r"The **dot plot** shows \1, with a panel each for \2",
    ),
    (
        re.compile(
            re.escape("The four **dot plots** show the selected genes in each cell subgroup for R and "
                      "NR patients, before treatment (top) and after treatment (bottom). Compare left "
                      "with right for R vs NR, and top with bottom for the effect of treatment.")
        ),
        "The **dot plot** shows the selected genes in each cell subgroup, with one panel per group: "
        "R and NR at baseline, then R and NR post-treatment. Compare neighbouring panels for R vs "
        "NR, and the first pair with the second for the effect of treatment.",
    ),
    (
        re.compile(
            re.escape("The **dot plots** show the selected genes in each cell subgroup at baseline "
                      "and post-treatment.")
        ),
        "The **dot plot** shows the selected genes in each cell subgroup, with a panel for baseline "
        "and one for post-treatment.",
    ),
    (
        re.compile(re.escape("The **dot plots** show the selected genes in each cell group.")),
        "The **dot plot** shows the selected genes in each cell group, with a panel each for "
        "responders and non-responders.",
    ),
    (
        re.compile(re.escape("The **dot plots** show the selected genes for each value of ")),
        "The **dot plot** shows the selected genes, with a panel for each value of ",
    ),
]

# Wording that means a text box still describes separate dot plots.
STALE_TEXT = ("**dot plots** show", "dot plot's disease")


def split_column(chart: dict[str, Any]) -> str | None:
    """The background filter column of a legacy dot plot restricted to some values, else None."""
    if chart.get("type") != LEGACY:
        return None
    bf = chart.get("background_filter")
    if not isinstance(bf, dict):
        return None
    column = bf.get("column")
    categories = bf.get("categories")
    if not isinstance(column, str) or column in ("", NO_COLUMN):
        return None
    if not isinstance(categories, list) or not categories:
        return None
    return column


def find_groups(charts: list[dict[str, Any]]) -> list[list[int]]:
    """Indexes of legacy dot plots that differ only in their background filter values."""
    groups: dict[tuple[str, str], list[int]] = {}
    for index, chart in enumerate(charts):
        column = split_column(chart)
        if column is not None:
            key = (json.dumps(chart.get("param"), sort_keys=True), column)
            groups.setdefault(key, []).append(index)
    return [indexes for indexes in groups.values() if len(indexes) > 1]


def join_words(words: list[str]) -> str:
    return words[0] if len(words) == 1 else ", ".join(words[:-1]) + " and " + words[-1]


def merged_title(titles: list[str], column_label: str) -> str:
    """Keep the shared start of the titles, e.g. "Selected genes in RA and PsA"."""
    prefix = os.path.commonprefix(titles)
    prefix = prefix[: prefix.rfind(" ") + 1]
    if prefix.strip() and all(len(title) > len(prefix) for title in titles):
        return prefix + join_words([title[len(prefix):] for title in titles])
    return f"Selected genes by {column_label}"


def box(chart: dict[str, Any]) -> tuple[float, float, float, float] | None:
    position, size = chart.get("position"), chart.get("size")
    if not (isinstance(position, list) and isinstance(size, list) and len(position) == 2 and len(size) == 2):
        return None
    x, y = position
    w, h = size
    if not all(isinstance(v, (int, float)) for v in (x, y, w, h)):
        return None
    return x, y, x + w, y + h


def merge(charts: list[dict[str, Any]], column_label: str) -> dict[str, Any]:
    """One compound dot plot covering `charts`, with their filter values as panels."""
    boxes = [b for b in map(box, charts) if b is not None]
    if len(boxes) == len(charts):
        # panels read left to right, top to bottom, as the separate plots were laid out
        charts = [c for _, c in sorted(zip(boxes, charts), key=lambda pair: (pair[0][1], pair[0][0]))]
    first = charts[0]
    values: list[str] = []
    for chart in charts:
        for value in chart["background_filter"]["categories"]:
            if value not in values:
                values.append(value)

    merged: dict[str, Any] = {key: first[key] for key in KEPT_KEYS if key in first}
    merged["type"] = COMPOUND
    merged["title"] = merged_title([str(c.get("title", "")) for c in charts], column_label)
    merged["split_by"] = first["background_filter"]["column"]
    merged["split_values"] = values
    old_scale = first.get("color_scale")
    scale: dict[str, Any] = old_scale if isinstance(old_scale, dict) else {}
    palette = scale.get("palette") if scale.get("palette") in PALETTES else PALETTES[0]
    merged["color_scale"] = {"palette": palette, "log": bool(scale.get("log", False)), "mode": "auto"}
    for legend in ("color_legend", "fraction_legend"):
        old = first.get(legend)
        merged[legend] = {"display": old.get("display", True) if isinstance(old, dict) else True}
    merged["reserve_legend_space"] = True
    axis = first.get("axis")
    if isinstance(axis, dict):
        merged["axis"] = {name: axis[name] for name in ("x", "y") if isinstance(axis.get(name), dict)}
    if boxes:
        left, top = min(b[0] for b in boxes), min(b[1] for b in boxes)
        merged["position"] = [left, top]
        merged["size"] = [max(b[2] for b in boxes) - left, max(b[3] for b in boxes) - top]
    for key in ("gsposition", "gssize"):
        if key in first:
            merged[key] = first[key]
    return merged


def update_text(text: str) -> tuple[str, int]:
    count = 0
    for pattern, replacement in TEXT_UPDATES:
        text, n = pattern.subn(replacement, text)
        count += n
    return text, count


def column_labels(datasources: list[dict[str, Any]]) -> dict[tuple[str, str], str]:
    labels: dict[tuple[str, str], str] = {}
    for ds in datasources:
        for column in ds.get("columns", []):
            field = column.get("field")
            if isinstance(field, str):
                labels[(ds.get("name", ""), field)] = str(column.get("name") or field)
    return labels


def convert_views(views: dict[str, Any], labels: dict[tuple[str, str], str]) -> list[dict[str, Any]]:
    """Merge dot plot groups in place; returns what changed, per view."""
    changes: list[dict[str, Any]] = []
    for view_name, view in views.items():
        charts_by_ds = view.get("initialCharts") if isinstance(view, dict) else None
        if not isinstance(charts_by_ds, dict):
            continue
        merges: list[dict[str, Any]] = []
        for ds_name, charts in charts_by_ds.items():
            if not isinstance(charts, list):
                continue
            groups = find_groups(charts)
            if not groups:
                continue
            replaced: dict[int, dict[str, Any]] = {}
            dropped: set[int] = set()
            for indexes in groups:
                column = str(charts[indexes[0]]["background_filter"]["column"])
                merged = merge([charts[i] for i in indexes], labels.get((ds_name, column), column))
                replaced[indexes[0]] = merged
                dropped.update(indexes[1:])
                merges.append({
                    "datasource": ds_name,
                    "merged": [str(charts[i].get("title", "")) for i in indexes],
                    "into": merged["title"],
                    "split_by": column,
                    "split_values": merged["split_values"],
                })
            charts_by_ds[ds_name] = [
                replaced.get(i, chart) for i, chart in enumerate(charts) if i not in dropped
            ]
        if not merges:
            continue
        text_updates = 0
        stale: list[str] = []
        for charts in charts_by_ds.values():
            for chart in charts if isinstance(charts, list) else []:
                if chart.get("type") == "text_box_chart" and isinstance(chart.get("text"), str):
                    chart["text"], n = update_text(chart["text"])
                    text_updates += n
                    if any(phrase in chart["text"] for phrase in STALE_TEXT):
                        stale.append(str(chart.get("title", "text box")))
        changes.append({"view": view_name, "merges": merges, "text_updates": text_updates, "stale_text": stale})
    return changes


def atomic_write_json(path: Path, data: Any) -> None:
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.", suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(data, handle, indent=2)
            handle.write("\n")
        os.replace(tmp, path)
    except BaseException:
        Path(tmp).unlink(missing_ok=True)
        raise


def convert_project(project: Path, *, write: bool = False) -> dict[str, Any]:
    views_path = project / "views.json"
    views = json.loads(views_path.read_text(encoding="utf-8"))
    datasources_path = project / "datasources.json"
    datasources = json.loads(datasources_path.read_text(encoding="utf-8")) if datasources_path.exists() else []
    changes = convert_views(views, column_labels(datasources))
    backup = None
    if write and changes:
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        backup = views_path.with_name(f"views.json.bak.{stamp}")
        backup.write_bytes(views_path.read_bytes())
        atomic_write_json(views_path, views)
    return {"project": str(project), "changes": changes, "written": bool(write and changes), "backup": str(backup) if backup else None}


def print_summary(summary: dict[str, Any]) -> None:
    print(f"Project: {summary['project']}")
    if not summary["changes"]:
        print("  Nothing to merge.")
    for change in summary["changes"]:
        print(f"  View: {change['view']}")
        for item in change["merges"]:
            print(f"    {len(item['merged'])} dot plots -> \"{item['into']}\" "
                  f"split by {item['split_by']}: {', '.join(item['split_values'])}")
        print(f"    Text updates: {change['text_updates']}")
        for title in change["stale_text"]:
            print(f"    Check by hand: \"{title}\" still describes separate dot plots")
    if summary["written"]:
        print(f"  Written. Backup: {summary['backup']}")
    elif summary["changes"]:
        print("  Dry run: nothing written. Add --write to apply.")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--project", type=Path, required=True, help="MDV project folder (has views.json)")
    parser.add_argument("--write", action="store_true", help="Apply the changes (default: dry run)")
    args = parser.parse_args(argv)
    try:
        summary = convert_project(args.project.expanduser(), write=args.write)
    except (OSError, ValueError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    print_summary(summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
