import pytest

from mdvtools.jobs.registry import CONCAT_COLUMNS, get_tool, validate_params


class FakeProject:
    """Minimal stand-in: validate_params only calls get_datasource_metadata"""

    def __init__(self, datasources):
        self._ds = datasources  # {ds_name: [field, ...]}

    def get_datasource_metadata(self, name):
        return {"columns": [{"field": f} for f in self._ds[name]]}


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
