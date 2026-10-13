from dataclasses import dataclass, asdict
from typing import Any


@dataclass(frozen=True)
class ParamSpec:
    name: str
    type: str  # "datasource" | "subgroup" | "column" | "selection" | "dropdown" | "text" | "int" | "float"
    label: str
    options_from: str | None = None  # column or subgroup picker scoped to THIS datasource param
    default: Any = None
    applies_to: str | None = None # compute param -> which worker call it feeds; None = a control/output param (eg. output_name)

@dataclass(frozen=True)
class OutputSpec:
    shape: str  # "column" -> column(ds, cols)
    datasource_param: str
    columns_param: str


@dataclass(frozen=True)
class ToolSpec:
    id: str
    name: str
    description: str
    params: list[ParamSpec]
    output: OutputSpec
    entrypoint: str  # "module: fn" - registry doubles as a dispatch table
    input_shape: str  # "columns" | "matrix"


CONCAT_COLUMNS = ToolSpec(
    id="concat_columns",
    name="Concatenate Columns",
    description="Concatenate multiple columns into a single column.",
    params=[
        ParamSpec("datasource", "datasource", "Datasource"),
        ParamSpec("column_a", "column", "First column", options_from="datasource"),
        ParamSpec("column_b", "column", "Second column", options_from="datasource"),
        ParamSpec("separator", "text", "Separator", default="_"),
        ParamSpec("output_name", "text", "New column name"),
    ],
    output=OutputSpec("column", "datasource", "output_name"),
    entrypoint="mdvtools.jobs.workers.concat_worker:run",
    input_shape="columns",
)

UMAP = ToolSpec(
    id="umap",
    name="UMAP",
    description="Compute a UMAP embedding from the expression matrix (adds UMAP_1, UMAP_2, ..., UMAP_n)",
    params=[
        ParamSpec("datasource", "datasource", "Datasource"),
        ParamSpec("layer", "subgroup", "Matrix", options_from="datasource", default="gs"),
        # the selection's name; empty runs on all rows (ADR-0017)
        ParamSpec("selection", "selection", "Cells", options_from="datasource"),
        ParamSpec("output_name", "text", "New column base name", default="UMAP"),
        ParamSpec("n_comps", "int", "PCA components", default=50, applies_to="pca"),
        ParamSpec("n_neighbors", "int", "Neighbors", default=15, applies_to="neighbors"),
        ParamSpec("min_dist", "float", "Minimum distance", default=0.5, applies_to="umap"),
        ParamSpec("n_components", "int", "Dimensions", default=2, applies_to="umap"),
    ],
    output=OutputSpec("column", "datasource", "output_name"),
    entrypoint="mdvtools.jobs.workers.umap_worker:run",
    input_shape="matrix",
)

REGISTRY: dict[str, ToolSpec] = {CONCAT_COLUMNS.id: CONCAT_COLUMNS, UMAP.id: UMAP}

def serialize_registry() -> list[dict]:
    """
    Client facing view of the tool registry for GET /jobs/tools.

    asdict each ToolSpec into the JSON the selector renders
    """
    tools = []
    for spec in REGISTRY.values():
        d = asdict(spec)
        d.pop("entrypoint", None)       # client doesn't need entrypoint
        tools.append(d)
    return tools


def get_tool(tool_id: str) -> ToolSpec:
    if tool_id not in REGISTRY:
        raise KeyError(f"Unknown tool: {tool_id}")
    return REGISTRY[tool_id]


def validate_params(spec: ToolSpec, params: dict, project, rows: list[int] | None = None) -> None:
    """ADR: 0006 backend re-validates. A "datasource" param must name a datasource; every "column"
    or "subgroup" param must name a column or matrix of the datasource its options_from points at.
    `rows` are the indices of a selection (ADR-0017): only a tool with a "selection" param takes
    them, a named selection needs them, and each must be a distinct row of the datasource."""
    if rows is not None and not any(p.type == "selection" for p in spec.params):
        raise ValueError(f"{spec.id} cannot run on a selection")
    for p in spec.params:
        if p.type == "datasource":
            if params.get(p.name) not in project.get_datasource_names():
                raise ValueError(f"{params.get(p.name)!r} is not a datasource")
        elif p.type == "subgroup":
            ds_name = params[p.options_from]
            keys = {
                k
                for ln in project.get_links(ds_name, "rows_as_columns")
                for k in (ln["link"].get("rows_as_columns") or {}).get("subgroups") or {}
            }
            if params.get(p.name) not in keys:
                raise ValueError(f"{params.get(p.name)!r} is not a matrix of {ds_name!r}")
        elif p.type == "column":
            ds_name = params[p.options_from]
            fields = {
                c["field"] for c in project.get_datasource_metadata(ds_name)["columns"]
            }
            if params.get(p.name) not in fields:
                raise ValueError(
                    f"{params.get(p.name)!r} is not a column of {ds_name!r}"
                )
        elif p.type == "selection":
            ds_name = params[p.options_from]
            if not params.get(p.name):
                if rows is not None:
                    raise ValueError("a selection needs a name")
            elif rows is None:
                raise ValueError(f"selection {params[p.name]!r} has no rows")
            elif not rows:
                raise ValueError("selection rows must not be empty")
            elif not all(isinstance(r, int) and not isinstance(r, bool) for r in rows):
                raise ValueError("selection rows must be whole numbers")
            else:
                size = int(project.get_datasource_metadata(ds_name)["size"])
                if min(rows) < 0 or max(rows) >= size:
                    raise ValueError(f"selection rows fall outside {ds_name!r}, which has {size} rows")
                if len(set(rows)) != len(rows):
                    raise ValueError("selection rows contain duplicates")
        elif p.type in ("int", "float") and params.get(p.name) is not None:
            v = params[p.name]
            ok = (not isinstance(v, bool)) and (
                isinstance(v, int) if p.type == "int" else isinstance(v, (int, float))
            )
            if not ok:
                raise ValueError(f"{p.name} must be {p.type}, got {type(v).__name__}")

    if not params.get(spec.output.columns_param):
        raise ValueError(f"{spec.output.columns_param} is required.")
