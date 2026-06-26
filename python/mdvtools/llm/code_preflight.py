"""
Preflight validation for ChatMDV-generated Python code.

This module performs AST-driven checks before subprocess execution to catch
deterministic failures (for example, missing chart imports) early.
"""
from __future__ import annotations

import ast
import importlib
import inspect
import re
from dataclasses import dataclass
from functools import lru_cache
from typing import Any

from mdvtools.llm.column_field_resolve import expression_wrapper_subgroup_key

_WRAPPER_RE = re.compile(r"([^|]+)\|([^|(]+)\(\1\)\|\s*(\d+)\s*$")

_FORBIDDEN_EXPR_ATTRS = frozenset({"feature_table", "feature_datasource"})

# Known hallucinated chart APIs (not defined on any mdvtools.charts class).
_FORBIDDEN_CHART_METHODS = frozenset({"set_row_indices", "set_background_filter"})

_NAMES_REQUIRING_IMPORT = frozenset(
    {
        "infer_datasource_roles",
        "build_expression_wrapper_token",
        "categorical_field_ids_from_metadata",
    }
)

# Pandas/AnnData-style key; MDV metadata column dicts use ``datatype``.
_FORBIDDEN_METADATA_COLUMN_KEYS = frozenset({"dtype"})


CHART_CLASS_MODULES: dict[str, str] = {
    "DotPlot": "mdvtools.charts.dot_plot",
    "ScatterPlot": "mdvtools.charts.scatter_plot",
    "HeatmapPlot": "mdvtools.charts.heatmap_plot",
    "HistogramPlot": "mdvtools.charts.histogram_plot",
    "BoxPlot": "mdvtools.charts.box_plot",
    "ViolinPlot": "mdvtools.charts.violin_plot",
    "DensityScatterPlot": "mdvtools.charts.density_scatter_plot",
    "ScatterPlot3D": "mdvtools.charts.scatter_plot_3D",
    "RowChart": "mdvtools.charts.row_chart",
    "StackedRowChart": "mdvtools.charts.stacked_row_plot",
    "PieChart": "mdvtools.charts.pie_chart",
    "RingChart": "mdvtools.charts.ring_chart",
    "AbundanceBoxPlot": "mdvtools.charts.abundance_box_plot",
    "MultiLinePlot": "mdvtools.charts.multi_line_plot",
    "TablePlot": "mdvtools.charts.table_plot",
    "WordcloudPlot": "mdvtools.charts.wordcloud_plot",
    "SankeyPlot": "mdvtools.charts.sankey_plot",
    "SelectionDialogPlot": "mdvtools.charts.selection_dialog_plot",
    "RowSummaryBox": "mdvtools.charts.row_summary_box_plot",
    "TextBox": "mdvtools.charts.text_box_plot",
}

_CHART_LIKE_SUFFIXES = ("Chart", "Plot", "Box")
_CHART_NAME_HINTS: dict[str, str] = {
    "Histogram": "HistogramPlot",
    "MultilineChart": "MultiLinePlot",
    "MultiLineChart": "MultiLinePlot",
}


@dataclass(frozen=True)
class PreflightIssue:
    code: str
    message: str
    line: int | None = None
    column: int | None = None
    symbol: str | None = None


@dataclass(frozen=True)
class PreflightResult:
    ok: bool
    issues: list[PreflightIssue]
    chart_classes: list[str]

    @property
    def summary(self) -> str:
        if not self.issues:
            return "Preflight passed."
        return "; ".join(issue.message for issue in self.issues)


def _collect_defined_names(tree: ast.AST) -> set[str]:
    defined: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                defined.add(alias.asname or alias.name.split(".")[0])
        elif isinstance(node, ast.ImportFrom):
            for alias in node.names:
                if alias.name == "*":
                    continue
                defined.add(alias.asname or alias.name)
        elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            defined.add(node.name)
        elif isinstance(node, ast.Assign):
            for t in node.targets:
                if isinstance(t, ast.Name):
                    defined.add(t.id)
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            defined.add(node.target.id)
    return defined


def _collect_imported_modules(tree: ast.AST) -> dict[str, str]:
    imports: dict[str, str] = {}
    for node in ast.walk(tree):
        if not isinstance(node, ast.ImportFrom):
            continue
        if not node.module:
            continue
        for alias in node.names:
            if alias.name == "*":
                continue
            imports[alias.asname or alias.name] = node.module
    return imports


def _extract_chart_calls(tree: ast.AST) -> list[tuple[str, ast.Call]]:
    calls: list[tuple[str, ast.Call]] = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        fn = node.func
        if isinstance(fn, ast.Name):
            calls.append((fn.id, node))
        elif isinstance(fn, ast.Attribute):
            calls.append((fn.attr, node))
    return calls


def _is_chart_like_name(name: str) -> bool:
    if name in CHART_CLASS_MODULES:
        return True
    if name in _CHART_NAME_HINTS:
        return True
    if not name or not name[0].isupper():
        return False
    return name.endswith(_CHART_LIKE_SUFFIXES)


def _build_unknown_chart_message(name: str) -> str:
    hint = _CHART_NAME_HINTS.get(name)
    if hint is not None:
        return (
            f"Unknown chart class `{name}`. "
            f"Use an allowlisted mdvtools chart class (did you mean `{hint}`?)."
        )
    return (
        f"Unknown chart class `{name}`. "
        "Use an allowlisted mdvtools chart class."
    )


def _extract_string_list(node: ast.AST) -> list[str] | None:
    if not isinstance(node, (ast.List, ast.Tuple)):
        return None
    values: list[str] = []
    for elt in node.elts:
        if isinstance(elt, ast.Constant) and isinstance(elt.value, str):
            values.append(elt.value)
        else:
            return None
    return values


def _collect_literal_bindings(tree: ast.AST) -> dict[str, str | list[str]]:
    bindings: dict[str, str | list[str]] = {}
    for node in ast.walk(tree):
        if not isinstance(node, ast.Assign):
            continue
        if len(node.targets) != 1 or not isinstance(node.targets[0], ast.Name):
            continue
        name = node.targets[0].id
        value = node.value
        if isinstance(value, ast.Constant) and isinstance(value.value, str):
            bindings[name] = value.value
            continue
        maybe_list = _extract_string_list(value)
        if maybe_list is not None:
            bindings[name] = maybe_list
    return bindings


def _resolve_string_arg(node: ast.AST, bindings: dict[str, str | list[str]]) -> str | None:
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        return node.value
    if isinstance(node, ast.Name):
        v = bindings.get(node.id)
        if isinstance(v, str):
            return v
    return None


def _resolve_string_list_arg(
    node: ast.AST, bindings: dict[str, str | list[str]]
) -> list[str] | None:
    literal = _extract_string_list(node)
    if literal is not None:
        return literal
    if isinstance(node, ast.Name):
        v = bindings.get(node.id)
        if isinstance(v, list):
            return [str(x) for x in v]
    return None


def _extract_dataframe_column_requests(
    tree: ast.AST, bindings: dict[str, str | list[str]]
) -> list[tuple[str, list[str], ast.Call]]:
    requests: list[tuple[str, list[str], ast.Call]] = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        if not isinstance(node.func, ast.Attribute):
            continue
        if node.func.attr != "get_datasource_as_dataframe":
            continue

        datasource_name: str | None = None
        columns: list[str] | None = None

        if node.args:
            datasource_name = _resolve_string_arg(node.args[0], bindings)
        for kw in node.keywords:
            if kw.arg == "datasource":
                datasource_name = _resolve_string_arg(kw.value, bindings)
            elif kw.arg == "columns":
                columns = _resolve_string_list_arg(kw.value, bindings)

        if columns is None:
            continue
        if datasource_name is None:
            continue
        requests.append((datasource_name, columns, node))
    return requests


@lru_cache(maxsize=64)
def _load_constructor_signature(class_name: str) -> inspect.Signature | None:
    module_name = CHART_CLASS_MODULES.get(class_name)
    if not module_name:
        return None
    module = importlib.import_module(module_name)
    cls = getattr(module, class_name)
    return inspect.signature(cls.__init__)


def _check_keyword_args(class_name: str, call: ast.Call) -> list[PreflightIssue]:
    try:
        sig = _load_constructor_signature(class_name)
    except Exception as exc:  # pragma: no cover - defensive runtime safety
        return [
            PreflightIssue(
                code="constructor_lookup_failed",
                message=f"Could not inspect constructor for `{class_name}`: {exc}",
                line=getattr(call, "lineno", None),
                column=getattr(call, "col_offset", None),
                symbol=class_name,
            )
        ]
    if sig is None:
        return []

    has_var_kw = any(p.kind == inspect.Parameter.VAR_KEYWORD for p in sig.parameters.values())
    if has_var_kw:
        return []

    valid_kw = {
        name
        for name, param in sig.parameters.items()
        if name != "self"
        and param.kind in (inspect.Parameter.POSITIONAL_OR_KEYWORD, inspect.Parameter.KEYWORD_ONLY)
    }

    issues: list[PreflightIssue] = []
    for kw in call.keywords:
        if kw.arg is None:
            # **kwargs spread; skip static validation for this call.
            return []
        if kw.arg not in valid_kw:
            issues.append(
                PreflightIssue(
                    code="invalid_constructor_kwarg",
                    message=f"Invalid constructor kwarg `{kw.arg}` for `{class_name}`.",
                    line=getattr(call, "lineno", None),
                    column=getattr(call, "col_offset", None),
                    symbol=class_name,
                )
            )
    return issues


def _check_wrapper_and_roles_api(tree: ast.AST, defined_names: set[str]) -> list[PreflightIssue]:
    issues: list[PreflightIssue] = []
    imported = _collect_defined_names(tree)

    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            func = node.func
            if (
                isinstance(func, ast.Attribute)
                and func.attr == "get_datasource_roles"
                and isinstance(func.value, ast.Name)
                and func.value.id == "project"
            ):
                issues.append(
                    PreflightIssue(
                        code="hallucinated_project_api",
                        message=(
                            "`project.get_datasource_roles()` does not exist. "
                            "Use injected `CHATMDV_*` constants or `infer_datasource_roles(project)`."
                        ),
                        line=getattr(node, "lineno", None),
                        column=getattr(node, "col_offset", None),
                        symbol="get_datasource_roles",
                    )
                )

        if isinstance(node, ast.Attribute) and node.attr in _FORBIDDEN_EXPR_ATTRS:
            issues.append(
                PreflightIssue(
                    code="invalid_expression_role_attr",
                    message=(
                        f"`{node.attr}` is not a field on RowsAsColumnsExpression; "
                        "use `datasource_name`, `name_column`, or `CHATMDV_EXPR_DATASOURCE`."
                    ),
                    line=getattr(node, "lineno", None),
                    column=getattr(node, "col_offset", None),
                    symbol=node.attr,
                )
            )

        if isinstance(node, ast.Name) and node.id in _NAMES_REQUIRING_IMPORT:
            if node.id not in imported and node.id not in defined_names:
                ctx = getattr(node, "ctx", None)
                if isinstance(ctx, ast.Load):
                    issues.append(
                        PreflightIssue(
                            code="missing_roles_import",
                            message=(
                                f"`{node.id}` is used but not imported. "
                                "It is available from mdvtools.llm.datasource_roles / column_field_resolve "
                                "(also in the script preamble)."
                            ),
                            line=getattr(node, "lineno", None),
                            column=getattr(node, "col_offset", None),
                            symbol=node.id,
                        )
                    )

    return issues


def _collect_chart_instance_types(tree: ast.AST) -> dict[str, str]:
    """Map variable names to chart class names for ``var = ChartClass(...)`` assignments."""
    bindings: dict[str, str] = {}
    for node in ast.walk(tree):
        if not isinstance(node, ast.Assign):
            continue
        if len(node.targets) != 1 or not isinstance(node.targets[0], ast.Name):
            continue
        var_name = node.targets[0].id
        value = node.value
        if not isinstance(value, ast.Call):
            continue
        fn = value.func
        class_name: str | None = None
        if isinstance(fn, ast.Name) and fn.id in CHART_CLASS_MODULES:
            class_name = fn.id
        elif isinstance(fn, ast.Attribute) and fn.attr in CHART_CLASS_MODULES:
            class_name = fn.attr
        if class_name:
            bindings[var_name] = class_name
    return bindings


@lru_cache(maxsize=64)
def _load_chart_class(class_name: str) -> type[Any] | None:
    module_name = CHART_CLASS_MODULES.get(class_name)
    if not module_name:
        return None
    module = importlib.import_module(module_name)
    return getattr(module, class_name, None)


@lru_cache(maxsize=64)
def _allowed_chart_methods(class_name: str) -> frozenset[str]:
    cls = _load_chart_class(class_name)
    if cls is None:
        return frozenset()
    from mdvtools.charts.base_plot import BasePlot

    allowed: set[str] = set()
    for chart_cls in (BasePlot, cls):
        for name, member in inspect.getmembers(chart_cls, predicate=inspect.isfunction):
            if name.startswith("_"):
                continue
            allowed.add(name)
    return frozenset(allowed)


def _metadata_column_key_from_subscript(node: ast.Subscript) -> str | None:
    sl = node.slice
    if isinstance(sl, ast.Constant) and isinstance(sl.value, str):
        return sl.value
    return None


def _is_dtype_metadata_access(node: ast.AST) -> bool:
    """True when code reads ``dtype`` from a column metadata dict (not pandas ``.dtype``)."""
    if isinstance(node, ast.Subscript):
        key = _metadata_column_key_from_subscript(node)
        if key in _FORBIDDEN_METADATA_COLUMN_KEYS:
            return True
    if isinstance(node, ast.Call):
        func = node.func
        if (
            isinstance(func, ast.Attribute)
            and func.attr == "get"
            and len(node.args) >= 1
        ):
            arg0 = node.args[0]
            if isinstance(arg0, ast.Constant) and isinstance(arg0.value, str):
                if arg0.value in _FORBIDDEN_METADATA_COLUMN_KEYS:
                    return True
        for kw in node.keywords:
            if kw.arg is None:
                continue
            if kw.arg != "dtype":
                continue
            if isinstance(kw.value, ast.Constant) and kw.value.value == "dtype":
                return True
    return False


def _check_invalid_metadata_column_dtype_access(tree: ast.AST) -> list[PreflightIssue]:
    issues: list[PreflightIssue] = []
    for node in ast.walk(tree):
        if not _is_dtype_metadata_access(node):
            continue
        issues.append(
            PreflightIssue(
                code="invalid_metadata_column_key",
                message=(
                    "MDV datasource column metadata uses 'datatype', not 'dtype'. "
                    "Use CHATMDV_CATEGORICAL_FIELD_IDS or "
                    "categorical_field_ids_from_metadata("
                    "project.get_datasource_metadata(CHATMDV_OBS_DATASOURCE))."
                ),
                line=getattr(node, "lineno", None),
                column=getattr(node, "col_offset", None),
                symbol="dtype",
            )
        )
    return issues


def _check_chart_method_calls(
    tree: ast.AST,
    chart_bindings: dict[str, str],
) -> list[PreflightIssue]:
    issues: list[PreflightIssue] = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        func = node.func
        if not isinstance(func, ast.Attribute):
            continue
        method = func.attr
        if method in _FORBIDDEN_CHART_METHODS:
            issues.append(
                PreflightIssue(
                    code="forbidden_chart_method",
                    message=(
                        f"Chart method `{method}` does not exist in mdvtools.charts. "
                        "Do not subset charts by row index or invented filter setters; use "
                        "SelectionDialogPlot, a filtered datasource, or dataframe analysis "
                        "with print(...to_markdown()) for subset statistics."
                    ),
                    line=getattr(node, "lineno", None),
                    column=getattr(node, "col_offset", None),
                    symbol=method,
                )
            )
            continue
        if not isinstance(func.value, ast.Name):
            continue
        receiver = func.value.id
        class_name = chart_bindings.get(receiver)
        if not class_name:
            continue
        allowed = _allowed_chart_methods(class_name)
        if method in allowed:
            continue
        issues.append(
            PreflightIssue(
                code="unknown_chart_method",
                message=(
                    f"`{receiver}.{method}()` is not a method on `{class_name}` "
                    f"(or BasePlot). Use only documented chart setters such as "
                    f"{sorted(allowed)[:12]}{'...' if len(allowed) > 12 else ''}."
                ),
                line=getattr(node, "lineno", None),
                column=getattr(node, "col_offset", None),
                symbol=method,
            )
        )
    return issues


def _check_wrapper_subgroup_literals(
    tree: ast.AST,
    allowed_subgroup_keys: set[str] | None,
) -> list[PreflightIssue]:
    if not allowed_subgroup_keys:
        return []
    issues: list[PreflightIssue] = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Constant) or not isinstance(node.value, str):
            continue
        token = node.value
        if not _WRAPPER_RE.fullmatch(token):
            continue
        sg = expression_wrapper_subgroup_key(token)
        if sg is not None and sg not in allowed_subgroup_keys:
            issues.append(
                PreflightIssue(
                    code="invalid_wrapper_subgroup",
                    message=(
                        f"Wrapper token uses subgroup `{sg}` which is not in project metadata "
                        f"(allowed: {sorted(allowed_subgroup_keys)}). "
                        "Use `CHATMDV_EXPR_SUBGROUP_KEY` or `build_expression_wrapper_token`."
                    ),
                    line=getattr(node, "lineno", None),
                    column=getattr(node, "col_offset", None),
                    symbol=sg,
                )
            )
    return issues


def validate_generated_code_preflight(
    code: str,
    *,
    datasource_fields: dict[str, set[str]] | None = None,
    allowed_wrapper_subgroup_keys: set[str] | None = None,
) -> PreflightResult:
    issues: list[PreflightIssue] = []
    try:
        tree = ast.parse(code)
    except SyntaxError as exc:
        return PreflightResult(
            ok=False,
            issues=[
                PreflightIssue(
                    code="syntax_error",
                    message=f"Syntax error during preflight: {exc.msg}",
                    line=exc.lineno,
                    column=exc.offset,
                )
            ],
            chart_classes=[],
        )

    defined_names = _collect_defined_names(tree)
    imported_modules = _collect_imported_modules(tree)
    chart_calls = _extract_chart_calls(tree)
    chart_bindings = _collect_chart_instance_types(tree)
    literal_bindings = _collect_literal_bindings(tree)
    dataframe_requests = _extract_dataframe_column_requests(tree, literal_bindings)
    seen: list[str] = []

    for class_name, call in chart_calls:
        if class_name not in CHART_CLASS_MODULES:
            continue
        if class_name not in seen:
            seen.append(class_name)

        if isinstance(call.func, ast.Name) and class_name not in defined_names:
            issues.append(
                PreflightIssue(
                    code="missing_import",
                    message=f"`{class_name}` is used but not imported/defined.",
                    line=getattr(call, "lineno", None),
                    column=getattr(call, "col_offset", None),
                    symbol=class_name,
                )
            )
            continue

        expected_module = CHART_CLASS_MODULES[class_name]
        actual_module = imported_modules.get(class_name)
        if actual_module and actual_module != expected_module:
            issues.append(
                PreflightIssue(
                    code="unexpected_import_module",
                    message=f"`{class_name}` should be imported from `{expected_module}`, found `{actual_module}`.",
                    line=getattr(call, "lineno", None),
                    column=getattr(call, "col_offset", None),
                    symbol=class_name,
                )
            )

        issues.extend(_check_keyword_args(class_name, call))

    # Catch unknown chart-like constructor names before runtime NameError.
    for name, call in chart_calls:
        if _is_chart_like_name(name) and name not in CHART_CLASS_MODULES:
            issues.append(
                PreflightIssue(
                    code="unknown_chart_class",
                    message=_build_unknown_chart_message(name),
                    line=getattr(call, "lineno", None),
                    column=getattr(call, "col_offset", None),
                    symbol=name,
                )
            )

    issues.extend(_check_wrapper_and_roles_api(tree, defined_names))
    issues.extend(_check_invalid_metadata_column_dtype_access(tree))
    issues.extend(_check_chart_method_calls(tree, chart_bindings))
    issues.extend(_check_wrapper_subgroup_literals(tree, allowed_wrapper_subgroup_keys))

    if datasource_fields:
        for datasource_name, columns, call in dataframe_requests:
            available = datasource_fields.get(datasource_name)
            if available is None:
                continue
            missing = [c for c in columns if c not in available]
            if not missing:
                continue
            available_preview = sorted(list(available))[:30]
            hint_parts: list[str] = []
            for col in missing:
                owners = sorted(
                    ds
                    for ds, flds in datasource_fields.items()
                    if ds != datasource_name and col in flds
                )
                if owners:
                    hint_parts.append(
                        f"Column `{col}` exists on datasource(s): {owners}"
                    )
            hint_suffix = ""
            if hint_parts:
                hint_suffix = " " + " ".join(hint_parts) + "."
            issues.append(
                PreflightIssue(
                    code="unknown_datasource_column",
                    message=(
                        f"Datasource `{datasource_name}` does not contain column(s) {missing}. "
                        f"Available fields include: {available_preview}.{hint_suffix}"
                    ),
                    line=getattr(call, "lineno", None),
                    column=getattr(call, "col_offset", None),
                    symbol=datasource_name,
                )
            )

    return PreflightResult(ok=len(issues) == 0, issues=issues, chart_classes=seen)


def format_preflight_issues(issues: list[PreflightIssue]) -> str:
    if not issues:
        return "No preflight issues."
    lines: list[str] = []
    for idx, issue in enumerate(issues, start=1):
        loc = ""
        if issue.line is not None:
            loc = f" (line {issue.line})"
        lines.append(f"{idx}. [{issue.code}] {issue.message}{loc}")
    return "\n".join(lines)

