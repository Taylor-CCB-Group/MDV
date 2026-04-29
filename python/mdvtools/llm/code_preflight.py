"""
Preflight validation for ChatMDV-generated Python code.

This module performs AST-driven checks before subprocess execution to catch
deterministic failures (for example, missing chart imports) early.
"""
from __future__ import annotations

import ast
import importlib
import inspect
from dataclasses import dataclass
from functools import lru_cache
from typing import Any


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


def validate_generated_code_preflight(
    code: str,
    *,
    datasource_fields: dict[str, set[str]] | None = None,
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

    # Catch common hallucinations like MultilineChart vs MultiLinePlot.
    for name, call in chart_calls:
        if name.endswith("Chart") and name not in CHART_CLASS_MODULES:
            issues.append(
                PreflightIssue(
                    code="unknown_chart_class",
                    message=(
                        f"Unknown chart class `{name}`. "
                        "Use an allowlisted mdvtools chart class (for multiline, use `MultiLinePlot`)."
                    ),
                    line=getattr(call, "lineno", None),
                    column=getattr(call, "col_offset", None),
                    symbol=name,
                )
            )

    if datasource_fields:
        for datasource_name, columns, call in dataframe_requests:
            available = datasource_fields.get(datasource_name)
            if available is None:
                continue
            missing = [c for c in columns if c not in available]
            if not missing:
                continue
            available_preview = sorted(list(available))[:30]
            issues.append(
                PreflightIssue(
                    code="unknown_datasource_column",
                    message=(
                        f"Datasource `{datasource_name}` does not contain column(s) {missing}. "
                        f"Available fields include: {available_preview}"
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

