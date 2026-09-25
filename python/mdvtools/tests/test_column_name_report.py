"""Tests for the standalone column-name inconsistency report."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

_SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "column_name_report.py"
_spec = importlib.util.spec_from_file_location("column_name_report", _SCRIPT)
if _spec is None or _spec.loader is None:
    raise RuntimeError(f"could not load {_SCRIPT}")
_report = importlib.util.module_from_spec(_spec)
sys.modules["column_name_report"] = _report
_spec.loader.exec_module(_report)

normalize_field = getattr(_report, "normalize_field")
discover_projects = getattr(_report, "discover_projects")
find_inconsistencies = getattr(_report, "find_inconsistencies")
scan_roots = getattr(_report, "scan_roots")
format_markdown = getattr(_report, "format_markdown")
write_csv = getattr(_report, "write_csv")


def _write_project(path: Path, datasources: list[dict[str, object]]) -> None:
    path.mkdir(parents=True, exist_ok=True)
    (path / "datasources.json").write_text(json.dumps(datasources), encoding="utf-8")
    (path / "nested").mkdir()
    (path / "nested" / "datasources.json").write_text(
        json.dumps([{"name": "cells", "columns": [{"field": "should_not_scan", "name": "should_not_scan"}]}]),
        encoding="utf-8",
    )


def test_normalize_joins_case_and_separators():
    assert normalize_field("cell_type") == "celltype"
    assert normalize_field("Cell Type") == "celltype"
    assert normalize_field("cell-type") == "celltype"
    assert normalize_field("celltype") == "celltype"


def test_identical_spellings_are_one_consistent_group():
    ColumnHit = getattr(_report, "ColumnHit")
    hits = [
        ColumnHit("/proj/a", "cells", "cell_type", "cell_type"),
        ColumnHit("/proj/b", "cells", "cell_type", "cell_type"),
    ]
    clusters = find_inconsistencies(hits)
    assert len(clusters) == 1
    assert clusters[0].inconsistent is False
    assert [spelling.field_name for spelling in clusters[0].spellings] == ["cell_type"]
    assert clusters[0].spellings[0].projects == ("/proj/a", "/proj/b")


def test_near_duplicate_fields_cluster_within_a_datasource(tmp_path: Path):
    _write_project(
        tmp_path / "a",
        [{"name": "cells", "columns": [{"field": "cell_type", "name": "Cell type"}]}],
    )
    _write_project(
        tmp_path / "b",
        [{"name": "cells", "columns": [{"field": "Cell Type", "name": "Cell Type"}]}],
    )
    _write_project(
        tmp_path / "c",
        [
            {
                "name": "cells",
                "columns": [
                    {"field": "celltype", "name": "celltype"},
                    {"field": "__index", "name": "__index"},
                ],
            },
            {"name": "genes", "columns": [{"field": "Cell Type", "name": "Cell Type"}]},
        ],
    )
    (tmp_path / "notes").mkdir()
    (tmp_path / "notes" / "readme.txt").write_text("not a project", encoding="utf-8")

    scan = scan_roots([tmp_path])
    assert len(scan.clusters) == 2
    cluster = next(item for item in scan.clusters if item.datasource == "cells")
    genes = next(item for item in scan.clusters if item.datasource == "genes")
    assert cluster.inconsistent is True
    assert cluster.normalized == ("celltype",)
    assert [spelling.field_name for spelling in cluster.spellings] == [
        "Cell Type",
        "cell_type",
        "celltype",
    ]
    assert genes.inconsistent is False
    assert [spelling.field_name for spelling in genes.spellings] == ["Cell Type"]
    fields_by_project = {
        spelling.field_name: spelling.projects for spelling in cluster.spellings
    }
    assert fields_by_project["cell_type"] == (str((tmp_path / "a").resolve()),)
    assert "should_not_scan" not in {
        spelling.field_name for item in scan.clusters for spelling in item.spellings
    }

    report = format_markdown(scan)
    assert "cell_type" in report
    assert "Cell Type" in report
    assert "celltype" in report
    assert "Cell type" in report
    assert "3 spellings" in report
    assert "1 spelling" in report
    assert "should_not_scan" not in report
    assert "__index" not in report


def test_fuzzy_merges_long_keys_and_skips_short_ones():
    ColumnHit = getattr(_report, "ColumnHit")
    hits = [
        ColumnHit("/proj/a", "cells", "cell_type", "cell_type"),
        ColumnHit("/proj/b", "cells", "cell_types", "cell_types"),
        ColumnHit("/proj/a", "cells", "pc1", "pc1"),
        ColumnHit("/proj/b", "cells", "pc2", "pc2"),
    ]
    exact = find_inconsistencies(hits)
    assert len(exact) == 4
    assert all(cluster.inconsistent is False for cluster in exact)
    fuzzy = find_inconsistencies(hits, fuzzy=0.8)
    assert len(fuzzy) == 3
    merged = next(cluster for cluster in fuzzy if len(cluster.spellings) > 1)
    assert merged.normalized == ("celltype", "celltypes")
    assert [spelling.field_name for spelling in merged.spellings] == [
        "cell_type",
        "cell_types",
    ]
    assert {cluster.normalized for cluster in fuzzy if not cluster.inconsistent} == {
        ("pc1",),
        ("pc2",),
    }


def test_columns_are_not_clustered_across_datasources():
    ColumnHit = getattr(_report, "ColumnHit")
    hits = [
        ColumnHit("/proj/a", "cells", "cell_type", "cell_type"),
        ColumnHit("/proj/b", "genes", "Cell Type", "Cell Type"),
    ]
    clusters = find_inconsistencies(hits)
    assert len(clusters) == 2
    assert {cluster.datasource for cluster in clusters} == {"cells", "genes"}
    assert all(cluster.inconsistent is False for cluster in clusters)


def test_load_omits_internal_fields_unless_requested(tmp_path: Path):
    _write_project(
        tmp_path / "a",
        [{"name": "cells", "columns": [{"field": "__index", "name": "__index"}]}],
    )
    _write_project(
        tmp_path / "b",
        [{"name": "cells", "columns": [{"field": "__Index", "name": "__Index"}]}],
    )
    hidden = scan_roots([tmp_path])
    assert hidden.clusters == []
    shown = scan_roots([tmp_path], include_internal=True)
    assert len(shown.clusters) == 1
    assert [spelling.field_name for spelling in shown.clusters[0].spellings] == [
        "__Index",
        "__index",
    ]


def test_discovery_ignores_non_projects_and_nested_datasources(tmp_path: Path):
    _write_project(
        tmp_path / "proj",
        [{"name": "cells", "columns": [{"field": "x", "name": "x"}]}],
    )
    (tmp_path / "loose").mkdir()
    (tmp_path / "loose" / "other.json").write_text("{}", encoding="utf-8")

    projects, errors = discover_projects([tmp_path])
    assert errors == []
    assert projects == [(tmp_path / "proj").resolve()]


def test_csv_rows_keep_display_name_with_its_project(tmp_path: Path):
    _write_project(
        tmp_path / "a",
        [{"name": "cells", "columns": [{"field": "cell_type", "name": "Cell type"}]}],
    )
    _write_project(
        tmp_path / "b",
        [{"name": "cells", "columns": [{"field": "celltype", "name": "celltype"}]}],
    )
    scan = scan_roots([tmp_path])
    destination = tmp_path / "report.csv"
    with destination.open("w", encoding="utf-8", newline="") as handle:
        write_csv(scan, handle)
    text = destination.read_text(encoding="utf-8")
    assert "datasource,normalized,field,display_name,project,inconsistent" in text
    assert "cell_type,Cell type," in text
    assert "celltype,celltype," in text
    assert ",true" in text
