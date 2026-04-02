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

if TYPE_CHECKING:
    from mdvtools.mdvproject import MDVProject

logger = logging.getLogger(__name__)

_WRAPPER_RE = re.compile(r"([^|]+)\|([^|(]+)\(\1\)\|\s*(\d+)")

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


def resolve_param_string(
    s: str,
    field_set: set[str],
    name_map: dict[str, str],
) -> str:
    """Map a single param string to canonical field id when needed."""
    if not isinstance(s, str):
        return s
    if _WRAPPER_RE.match(s):
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
    if "param" not in chart:
        return
    p = chart["param"]
    if isinstance(p, str):
        chart["param"] = resolve_param_string(p, field_set, name_map)
    elif isinstance(p, list):
        chart["param"] = [
            resolve_param_string(x, field_set, name_map) if isinstance(x, str) else x
            for x in p
        ]


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
    return out
