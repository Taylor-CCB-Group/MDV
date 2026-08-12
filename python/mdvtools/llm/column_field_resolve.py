"""
Resolve chart param strings from display name to datasource `field` id for Chat MDV.

The frontend indexes columns by `field`; LLM output may use `name`. Normalization
runs after generated code executes and updates the saved view JSON.
"""
from __future__ import annotations

import copy
import logging
from typing import TYPE_CHECKING, Any
import re

from mdvtools.llm.datasource_roles import CHAT_RANK_GENES_DATASOURCE_NAME

if TYPE_CHECKING:
    from mdvtools.mdvproject import MDVProject

logger = logging.getLogger(__name__)

_WRAPPER_RE = re.compile(r"([^|]+)\|([^|(]+)\(\1\)\|\s*(\d+)\s*$")
_MARKER_ALIASES = {
    "group": "cluster",
    "cluster": "group",
    "names": "gene",
    "gene": "names",
    "scores": "score",
    "score": "scores",
    "logfoldchanges": "logfoldchange",
    "logfoldchange": "logfoldchanges",
    "pvals": "pval",
    "pval": "pvals",
    "pvals_adj": "pval_adj",
    "pval_adj": "pvals_adj",
}

def build_name_to_field_map(columns: list[dict]) -> dict[str, str]:
    """
    Map display `name` -> `field` when unambiguous.
    If the same name maps to multiple distinct fields, that name is omitted.
    """
    name_to_fields: dict[str, set[str]] = {}
    for c in columns:
        name = c.get("name")
        if name is None:
            continue
        field = c.get("field", name)
        if field is None:
            continue
        name_to_fields.setdefault(str(name), set()).add(str(field))
    out: dict[str, str] = {}
    for name, fields in name_to_fields.items():
        if len(fields) == 1:
            out[name] = next(iter(fields))
        else:
            logger.debug(
                "Ambiguous column name %r -> multiple fields %s; skipping name map",
                name,
                fields,
            )
    return out


def field_set_from_columns(columns: list[dict]) -> set[str]:
    s: set[str] = set()
    for c in columns:
        f = c.get("field", c.get("name"))
        if f is not None:
            s.add(str(f))
    return s


def rows_as_columns_subgroup_keys_from_metadata(md: dict[str, Any]) -> set[str]:
    """
    Subgroup keys available for expression wrapper tokens on this datasource
    (from ``links.*.rows_as_columns.subgroups`` in datasource metadata).
    """
    keys: set[str] = set()
    links = md.get("links") or {}
    if not isinstance(links, dict):
        return keys
    for _target, link in links.items():
        if not isinstance(link, dict):
            continue
        rac = link.get("rows_as_columns") or {}
        if not isinstance(rac, dict):
            continue
        sgs = rac.get("subgroups") or {}
        if isinstance(sgs, dict):
            keys.update(str(k) for k in sgs.keys())
    return keys


def expression_wrapper_subgroup_key(token: str) -> str | None:
    """Return subgroup key from a wrapper token, or None if not a wrapper."""
    m = _WRAPPER_RE.fullmatch(token.strip()) if isinstance(token, str) else None
    if not m:
        return None
    return m.group(1).strip()


def build_expression_wrapper_token(subgroup_key: str, feature_label: str, index: int) -> str:
    """Canonical wrapper string for rows-as-columns expression columns (matches frontend/DataStore)."""
    sk = subgroup_key.strip()
    lab = feature_label.strip()
    return f"{sk}|{lab}({sk})|{int(index)}"


def resolve_param_string(
    s: str,
    field_set: set[str],
    name_map: dict[str, str],
) -> str:
    """Map a single param string to canonical field id when needed."""
    if not isinstance(s, str):
        return s
    if _WRAPPER_RE.fullmatch(s.strip()):
        return s
    if s in field_set:
        return s
    if s in name_map:
        return name_map[s]
    return s


def _normalize_chart_dict(
    chart: dict[str, Any],
    field_set: set[str],
    name_map: dict[str, str],
) -> None:
    if "param" in chart:
        p = chart["param"]
        if isinstance(p, str):
            chart["param"] = resolve_param_string(p, field_set, name_map)
        elif isinstance(p, list):
            chart["param"] = [
                resolve_param_string(x, field_set, name_map) if isinstance(x, str) else x
                for x in p
            ]
    _normalize_nested_column_refs(chart, field_set, name_map)


def _normalize_nested_column_refs(
    chart: dict[str, Any],
    field_set: set[str],
    name_map: dict[str, str],
) -> None:
    """Map display names to field ids on color_by / tooltip / background_filter when unambiguous."""
    cb = chart.get("color_by")
    if isinstance(cb, str):
        chart["color_by"] = resolve_param_string(cb, field_set, name_map)
    elif isinstance(cb, dict):
        col = cb.get("column")
        if isinstance(col, str):
            cb["column"] = resolve_param_string(col, field_set, name_map)
        elif isinstance(col, dict) and isinstance(col.get("field"), str):
            col["field"] = resolve_param_string(col["field"], field_set, name_map)

    tip = chart.get("tooltip")
    if isinstance(tip, dict):
        tc = tip.get("column")
        if isinstance(tc, str):
            tip["column"] = resolve_param_string(tc, field_set, name_map)
        elif isinstance(tc, list):
            tip["column"] = [
                resolve_param_string(x, field_set, name_map) if isinstance(x, str) else x
                for x in tc
            ]

    bf = chart.get("background_filter")
    if isinstance(bf, dict) and isinstance(bf.get("column"), str):
        bf["column"] = resolve_param_string(bf["column"], field_set, name_map)


def _normalize_marker_alias_tokens(
    chart: dict[str, Any], ds_name: str, field_set: set[str]
) -> None:
    """
    Normalize common Scanpy marker table aliases for chat marker datasource table charts.
    """
    if ds_name != CHAT_RANK_GENES_DATASOURCE_NAME:
        return
    if chart.get("type") != "table_chart":
        return
    p = chart.get("param")
    if not isinstance(p, list):
        return
    out: list[Any] = []
    for token in p:
        if not isinstance(token, str):
            out.append(token)
            continue
        if token in field_set:
            out.append(token)
            continue
        alt = _MARKER_ALIASES.get(token)
        if alt and alt in field_set:
            out.append(alt)
        else:
            out.append(token)
    chart["param"] = out


def _chart_param_string_tokens(chart: dict[str, Any]) -> list[str]:
    """Collect string tokens from chart ``param`` for validation."""
    p = chart.get("param")
    if p is None:
        return []
    if isinstance(p, str):
        return [p]
    if isinstance(p, list):
        return [str(x) for x in p if isinstance(x, str)]
    return []


def _chart_column_reference_strings(chart: dict[str, Any]) -> list[str]:
    """
    Collect datasource column field strings used by a chart config (all common chart types).

    Includes ``param`` plus legacy column slots like ``color_by``, ``tooltip.column``,
    ``background_filter.column`` when present as strings.
    """
    out: list[str] = []
    out.extend(_chart_param_string_tokens(chart))

    cb = chart.get("color_by")
    if isinstance(cb, str):
        out.append(cb)
    elif isinstance(cb, dict):
        col = cb.get("column")
        if isinstance(col, str):
            out.append(col)
        elif isinstance(col, dict):
            f = col.get("field")
            if isinstance(f, str):
                out.append(f)

    tip = chart.get("tooltip")
    if isinstance(tip, dict):
        tc = tip.get("column")
        if isinstance(tc, str):
            out.append(tc)
        elif isinstance(tc, list):
            out.extend(str(x) for x in tc if isinstance(x, str))

    bf = chart.get("background_filter")
    if isinstance(bf, dict):
        col = bf.get("column")
        if isinstance(col, str):
            out.append(col)

    return out


def _param_token_valid_for_datasource(
    token: str,
    field_set: set[str],
    subgroup_keys: set[str],
) -> bool:
    """
    True if ``token`` is a known column field, or a valid expression wrapper whose
    subgroup prefix exists in project metadata (no hallucinated subgroup keys).
    """
    if not isinstance(token, str):
        return False
    sg = expression_wrapper_subgroup_key(token)
    if sg is not None:
        if not subgroup_keys:
            return False
        return sg in subgroup_keys
    return token in field_set


def prune_view_charts_with_invalid_params(view: dict[str, Any], project: Any) -> tuple[dict[str, Any], int]:
    """
    Drop charts whose ``param`` strings are not valid field ids on that chart's datasource.

    Wrapper expression tokens (``subgroup|feature(subgroup)|index``) are kept only when
    ``subgroup`` matches a key under ``links.*.rows_as_columns.subgroups`` for that datasource.
    Used after :func:`normalize_view_chart_params` so display-name resolution has already run.
    """
    out = copy.deepcopy(view)
    initial = out.get("initialCharts")
    if not isinstance(initial, dict):
        return out, 0
    try:
        valid_ds = {str(d.get("name")) for d in (project.datasources or []) if d.get("name")}
    except Exception:
        valid_ds = set()

    removed = 0
    dropped_ds = 0
    for ds_name, charts in list(initial.items()):
        if valid_ds and ds_name not in valid_ds:
            dropped_ds += 1
            n_ch = len(charts) if isinstance(charts, list) else 0
            removed += n_ch
            logger.warning(
                "Removing initialCharts key %r: not a project datasource (have: %s)",
                ds_name,
                sorted(valid_ds) if valid_ds else "?",
            )
            del initial[ds_name]
            continue
        if not isinstance(charts, list):
            continue
        try:
            md = project.get_datasource_metadata(ds_name)
            columns = md.get("columns") or []
        except Exception as e:
            logger.warning("prune_view_charts_with_invalid_params: datasource %s: %s", ds_name, e)
            continue
        field_set = field_set_from_columns(columns)
        subgroup_keys = rows_as_columns_subgroup_keys_from_metadata(md)
        kept: list[Any] = []
        for ch in charts:
            if not isinstance(ch, dict):
                kept.append(ch)
                continue
            tokens = _chart_column_reference_strings(ch)
            bad = [
                t
                for t in tokens
                if not _param_token_valid_for_datasource(t, field_set, subgroup_keys)
            ]
            if bad:
                removed += 1
                logger.warning(
                    "Dropping chart type=%r on datasource %r: param tokens not in metadata: %s",
                    ch.get("type"),
                    ds_name,
                    bad,
                )
                continue
            kept.append(ch)
        initial[ds_name] = kept
    if dropped_ds:
        logger.warning(
            "prune_view_charts_with_invalid_params: dropped %d unknown datasource key(s)",
            dropped_ds,
        )
    return out, removed


def normalize_view_chart_params(
    view: dict[str, Any],
    project: Any,
) -> dict[str, Any]:
    """
    Return a deep copy of view with chart `param` strings resolved to field ids
    per datasource key in initialCharts (each datasource uses its own column metadata).
    """
    out = copy.deepcopy(view)
    initial = out.get("initialCharts")
    if not isinstance(initial, dict):
        return out

    for ds_name, charts in initial.items():
        if not isinstance(charts, list):
            continue
        try:
            md = project.get_datasource_metadata(ds_name)
            columns = md.get("columns") or []
        except Exception as e:
            logger.warning("normalize_view_chart_params: datasource %s: %s", ds_name, e)
            continue
        field_set = field_set_from_columns(columns)
        name_map = build_name_to_field_map(columns)
        for ch in charts:
            if isinstance(ch, dict):
                _normalize_chart_dict(ch, field_set, name_map)
                _normalize_marker_alias_tokens(ch, ds_name, field_set)
    return out


def ensure_view_gridstack_layout(view: dict[str, Any]) -> dict[str, Any]:
    """
    Ensure each datasource panel in a saved view uses gridstack layout.

    ChatMDV-generated views often only set ``initialCharts``; the frontend defaults
    missing layout to absolute positioning.
    """
    out = copy.deepcopy(view)
    initial = out.get("initialCharts")
    if not isinstance(initial, dict):
        return out

    data_sources = out.get("dataSources")
    if not isinstance(data_sources, dict):
        data_sources = {}
        out["dataSources"] = data_sources

    for ds_name in initial:
        panel = data_sources.get(ds_name)
        if not isinstance(panel, dict):
            panel = {}
            data_sources[ds_name] = panel
        panel["layout"] = "gridstack"

    return out
