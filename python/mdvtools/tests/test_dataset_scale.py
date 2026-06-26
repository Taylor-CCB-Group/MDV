import pytest

from mdvtools.llm.dataset_scale import (
    ProjectScale,
    agent_probe_columns_from_metadata,
    load_agent_dataframes,
)
from mdvtools.llm.datasource_roles import InferredDatasourceRoles, RowsAsColumnsExpression


class _ScaleProject:
    def __init__(self, *, metadata_by_name: dict[str, dict] | None = None):
        self._metadata_by_name = metadata_by_name or {}
        self.loads: list[tuple[str, list[str] | None]] = []

    def get_datasource_metadata(self, name: str) -> dict:
        return self._metadata_by_name.get(name, {"columns": [], "size": 200_000})

    def get_datasource_as_dataframe(self, name: str, columns: list[str] | None = None):
        self.loads.append((name, columns))
        return f"df:{name}"

    def get_datasource_names(self) -> list[str]:
        return ["table_a", "table_b", "table_c"]


_LARGE_SCALE = ProjectScale(
    obs_rows=200_000,
    obs_columns=10,
    estimated_obs_df_mb=30.0,
    available_ram_mb=8000.0,
    is_large=True,
    has_h5ad=False,
    obs_datasource="cells",
)

_OBS_COLUMNS = [{"field": "leiden", "datatype": "text"}]


def test_agent_probe_columns_returns_empty_when_metadata_has_no_columns():
    project = _ScaleProject()
    assert agent_probe_columns_from_metadata(project, "cells") == []


def test_load_agent_dataframes_large_obs_refuses_full_load_when_no_probe_columns():
    project = _ScaleProject()
    roles = InferredDatasourceRoles(obs_datasource="cells", expressions=[])

    with pytest.raises(ValueError, match="refusing full table load"):
        load_agent_dataframes(project, roles, _LARGE_SCALE)

    assert project.loads == []


def test_load_agent_dataframes_large_obs_uses_column_subset():
    project = _ScaleProject(metadata_by_name={"cells": {"columns": _OBS_COLUMNS}})
    roles = InferredDatasourceRoles(obs_datasource="cells", expressions=[])

    result = load_agent_dataframes(project, roles, _LARGE_SCALE)

    assert result == {"cells": "df:cells"}
    assert project.loads == [("cells", ["leiden"])]


def test_load_agent_dataframes_large_extra_skips_when_no_probe_columns():
    project = _ScaleProject(
        metadata_by_name={
            "cells": {"columns": _OBS_COLUMNS},
            "table_b": {"columns": []},
        }
    )
    roles = InferredDatasourceRoles(obs_datasource="cells", expressions=[])

    result = load_agent_dataframes(
        project, roles, _LARGE_SCALE, extra_datasource="table_b"
    )

    assert result == {"cells": "df:cells"}
    assert project.loads == [("cells", ["leiden"])]


def test_load_agent_dataframes_large_expr_skips_when_no_probe_columns():
    project = _ScaleProject(
        metadata_by_name={
            "cells": {"columns": _OBS_COLUMNS},
            "rna": {"columns": []},
        }
    )
    roles = InferredDatasourceRoles(
        obs_datasource="cells",
        expressions=[
            RowsAsColumnsExpression(
                datasource_name="rna",
                name_column="",
                subgroup_key="rna_expr",
                subgroup_label="rna_expr",
            )
        ],
    )

    result = load_agent_dataframes(project, roles, _LARGE_SCALE)

    assert result == {"cells": "df:cells"}
    assert project.loads == [("cells", ["leiden"])]


def test_load_agent_dataframes_selected_datasources_multi_table():
    project = _ScaleProject(
        metadata_by_name={
            "table_a": {"columns": _OBS_COLUMNS},
            "table_b": {"columns": [{"field": "assay", "datatype": "text"}]},
        }
    )
    roles = InferredDatasourceRoles(obs_datasource="table_a", expressions=[])

    result = load_agent_dataframes(
        project,
        roles,
        _LARGE_SCALE,
        selected_datasources=["table_a", "table_b"],
    )

    assert result == {"table_a": "df:table_a", "table_b": "df:table_b"}
    assert ("table_a", ["leiden"]) in project.loads
    assert ("table_b", ["assay"]) in project.loads
