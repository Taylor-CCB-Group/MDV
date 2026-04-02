from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Optional


@dataclass(frozen=True)
class RowsAsColumnsExpression:
    """
    A wrapper-capable expression datasource reachable via a rows-as-columns link.

    `subgroup_key` is the key used in the FieldName wrapper string:
        f"{subgroup_key}|{feature}({subgroup_key})|{row_index}"
    """

    datasource_name: str
    name_column: str
    subgroup_key: str
    subgroup_label: str


@dataclass(frozen=True)
class InferredDatasourceRoles:
    obs_datasource: str
    expressions: list[RowsAsColumnsExpression]

    def preferred_expression(self) -> Optional[RowsAsColumnsExpression]:
        if not self.expressions:
            return None
        # Prefer common modality names when present.
        for preferred in ("rna", "genes", "gene", "expression"):
            for e in self.expressions:
                if e.datasource_name.lower() == preferred:
                    return e
        return self.expressions[0]


def _pick_obs_datasource(project: Any) -> str:
    names = []
    try:
        names = list(project.get_datasource_names())
    except Exception:
        # Fallback: MDVProject.datasources property
        try:
            names = [ds.get("name") for ds in (project.datasources or []) if ds.get("name")]
        except Exception:
            names = []

    if "cells" in names:
        return "cells"
    # last-resort: first datasource
    if names:
        return str(names[0])
    raise ValueError("Project has no datasources; cannot infer roles")


def _pick_default_subgroup(subgroups: dict[str, Any]) -> tuple[str, str]:
    """
    Pick a reasonable default subgroup from rows-as-columns metadata.
    Prefer keys that look like primary expression (e.g. *_expr).
    """
    keys = list(subgroups.keys())
    if not keys:
        raise ValueError("rows_as_columns link has no subgroups")
    for k in keys:
        if k.lower().endswith("_expr") or k.lower() == "expr":
            sg = subgroups[k] or {}
            return k, str(sg.get("label") or sg.get("name") or k)
    k0 = keys[0]
    sg0 = subgroups[k0] or {}
    return k0, str(sg0.get("label") or sg0.get("name") or k0)


def infer_datasource_roles(project: Any) -> InferredDatasourceRoles:
    """
    Infer obs vs wrapper-expression datasources from MDV project metadata.

    This uses rows-as-columns links on the obs datasource to discover expression feature tables
    (e.g. `rna`, `protein`) and their available subgroups (e.g. `rna_expr`).
    """
    obs = _pick_obs_datasource(project)
    expressions: list[RowsAsColumnsExpression] = []

    # Discover rows-as-columns links from the obs datasource to other datasources.
    try:
        links = project.get_links(obs, "rows_as_columns")
    except Exception:
        links = []

    for item in links or []:
        try:
            target_ds = str(item.get("datasource"))
            link = item.get("link", {}).get("rows_as_columns", {})
            name_column = str(link.get("name_column") or "name")
            subgroups = link.get("subgroups") or {}
            if not isinstance(subgroups, dict) or len(subgroups) == 0:
                continue
            subgroup_key, subgroup_label = _pick_default_subgroup(subgroups)
            expressions.append(
                RowsAsColumnsExpression(
                    datasource_name=target_ds,
                    name_column=name_column,
                    subgroup_key=subgroup_key,
                    subgroup_label=subgroup_label,
                )
            )
        except Exception:
            # Best-effort inference; ignore malformed link entries.
            continue

    return InferredDatasourceRoles(obs_datasource=obs, expressions=expressions)

