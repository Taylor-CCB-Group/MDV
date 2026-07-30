import pytest

from mdvtools.jobs.registry import CONCAT_COLUMNS, get_tool, validate_params
from mdvtools.jobs.registry import ToolSpec, ParamSpec, OutputSpec

class FakeProject:
    """Minimal stand-in: validate_params only calls get_datasource_metadata"""

    def __init__(self, datasources):
        self._ds = datasources  # {ds_name: [field, ...]}

    def get_datasource_metadata(self, name):
        return {"columns": [{"field": f} for f in self._ds[name]]}


NUMERIC_SPEC = ToolSpec(
    id="numtool",
    name="Numeric tool",
    description="test",
    params=[
        ParamSpec("datasource", "dropdown", "Datasource"),
        ParamSpec("n_neighbors", "int", "Neighbors", default=15),
        ParamSpec("min_dist", "float", "Min distance", default=0.5),
        ParamSpec("output_name", "text", "Name", default="OUT"),
    ],
    output=OutputSpec("column", "datasource", "output_name"),
    entrypoint="x:run",
    input_shape="matrix",
)

@pytest.fixture
def project():
    return FakeProject({"cells": ["sample", "cluster"]})


def test_valid_params_pass(project):
    validate_params(
        CONCAT_COLUMNS,
        {
            "datasource": "cells",
            "column_a": "sample",
            "column_b": "cluster",
            "separator": "_",
            "output_name": "sample_cluster",
        },
        project,
    )  # no raise


def test_rejects_columns_not_on_datasource(project):
    with pytest.raises(ValueError, match="is not a column of 'cells'"):
        validate_params(
            CONCAT_COLUMNS,
            {
                "datasource": "cells",
                "column_a": "sample",
                "column_b": "nope",  # not a field of cells
                "output_name": "x",
            },
            project,
        )


def test_rejects_missing_output_name(project):
    with pytest.raises(ValueError, match="output_name is required"):
        validate_params(
            CONCAT_COLUMNS,
            {
                "datasource": "cells",
                "column_a": "sample",
                "column_b": "cluster",
                # output_name missing
            },
            project,
        )


def test_get_tool_unknown_raises():
    with pytest.raises(KeyError, match="Unknown tool: "):
        get_tool("unknown")


def test_numeric_params_absent_falls_back(project):
    # absent numeric params are fine — the worker uses scanpy's defaults
    validate_params(NUMERIC_SPEC, {"datasource": "cells", "output_name": "OUT"}, project)


def test_numeric_params_correct_type_pass(project):
    validate_params(
        NUMERIC_SPEC,
        {"datasource": "cells", "n_neighbors": 30, "min_dist": 0.1, "output_name": "OUT"},
        project,
    )


def test_float_param_accepts_int(project):
    # an int is a valid float
    validate_params(NUMERIC_SPEC, {"datasource": "cells", "min_dist": 1, "output_name": "OUT"}, project)


def test_rejects_non_int_for_int_param(project):
    with pytest.raises(ValueError, match="n_neighbors must be int"):
        validate_params(
            NUMERIC_SPEC, {"datasource": "cells", "n_neighbors": 3.5, "output_name": "OUT"}, project
        )


def test_rejects_string_for_float_param(project):
    with pytest.raises(ValueError, match="min_dist must be float"):
        validate_params(
            NUMERIC_SPEC, {"datasource": "cells", "min_dist": "0.1", "output_name": "OUT"}, project
        )


def test_rejects_bool_for_int_param(project):
    with pytest.raises(ValueError, match="n_neighbors must be int"):
        validate_params(
            NUMERIC_SPEC, {"datasource": "cells", "n_neighbors": True, "output_name": "OUT"}, project
        )
