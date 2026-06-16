from pathlib import Path

import pandas
import pytest

pytest.importorskip("pyranges")

from mdvtools.mdvproject import MDVProject
from mdvtools.feature_mapper import add_features_to_project, map_sv_genes


def _write_bed(path: Path, rows: list[tuple[str, int, int, str]]) -> str:
    path.write_text(
        "".join(f"{chrom}\t{start}\t{end}\t{name}\n" for chrom, start, end, name in rows),
        encoding="utf-8",
    )
    return str(path)


def test_map_sv_genes_breakpoint_and_range_modes(tmp_path):
    genes_bed = _write_bed(
        tmp_path / "genes.bed",
        [
            ("chr1", 100, 200, "GENE_A"),
            ("chr1", 300, 360, "GENE_B"),
            ("chr2", 395, 430, "GENE_C"),
        ],
    )
    sv_df = pandas.DataFrame(
        {
            "chr1": ["chr1", "chr1"],
            "pos1": [120, 320],
            "chr2": ["chr1", "chr2"],
            "pos2": [180, 410],
            "svtype": ["DEL", "BND"],
        },
        index=pandas.Index(["sv_del", "sv_bnd"]),
    )

    result = map_sv_genes(sv_df, genes_bed, window_bp=10).set_index("sv_id")

    assert result.loc["sv_del", "genes"] == "GENE_A"
    assert result.loc["sv_bnd", "genes"] == "GENE_B,GENE_C"


def test_add_features_to_project_interval_datasource(tmp_path):
    project = MDVProject(str(tmp_path / "proj"), delete_existing=True)
    interval_df = pandas.DataFrame(
        {
            "chromosome": ["1", "chr1", "chr2"],
            "left": [100, 180, 10],
            "right": [220, 260, 20],
        }
    )
    project.add_datasource("intervals", interval_df)
    metadata = project.get_datasource_metadata("intervals")
    metadata["genome"] = {
        "type": "interval",
        "assembly": "hg38",
        "columns": {"chr": "chromosome", "start": "left", "end": "right"},
    }
    project.set_datasource_metadata(metadata)

    genes_bed = _write_bed(
        tmp_path / "interval_features.bed",
        [
            ("chr1", 90, 160, "GENE_A"),
            ("chr1", 150, 240, "GENE_B"),
            ("chr2", 50, 100, "GENE_C"),
        ],
    )

    add_features_to_project(project, "intervals", genes_bed)

    assert project.get_column("intervals", "features") == ["GENE_A,GENE_B", "GENE_B", ""]
    column_meta = project.get_column_metadata("intervals", "features")
    assert column_meta["datatype"] == "multitext"
    assert column_meta["delimiter"] == ","


def test_add_features_to_project_sv_datasource(tmp_path):
    project = MDVProject(str(tmp_path / "proj"), delete_existing=True)
    sv_df = pandas.DataFrame(
        {
            "chrom_a": ["chr1", "chr1", "chr4"],
            "start_a": [105, 320, 1],
            "chrom_b": ["chr1", "chr2", "chr4"],
            "start_b": [205, 408, 1],
            "kind": ["DEL", "BND", "DEL"],
            "span": [100, 88, 0],
        }
    )
    project.add_datasource("svs", sv_df)
    metadata = project.get_datasource_metadata("svs")
    metadata["genome"] = {
        "type": "sv",
        "assembly": "hg38",
        "columns": {
            "chr1": "chrom_a",
            "pos1": "start_a",
            "chr2": "chrom_b",
            "pos2": "start_b",
            "svtype": "kind",
            "length": "span",
        },
    }
    project.set_datasource_metadata(metadata)

    genes_bed = _write_bed(
        tmp_path / "sv_features.bed",
        [
            ("chr1", 100, 210, "GENE_A"),
            ("chr1", 315, 330, "GENE_B"),
            ("chr2", 400, 420, "GENE_C"),
        ],
    )

    add_features_to_project(project, "svs", genes_bed, window_bp=5)

    assert project.get_column("svs", "features") == ["GENE_A", "GENE_B,GENE_C", ""]
    assert project.get_column_metadata("svs", "features")["datatype"] == "multitext"


def test_add_features_to_project_rejects_missing_or_invalid_genome_metadata(tmp_path):
    project = MDVProject(str(tmp_path / "proj"), delete_existing=True)
    project.add_datasource("ds", pandas.DataFrame({"chr": ["chr1"], "start": [1], "end": [2]}))

    with pytest.raises(ValueError, match="missing genome metadata"):
        add_features_to_project(project, "ds", str(tmp_path / "missing.bed"))

    metadata = project.get_datasource_metadata("ds")
    metadata["genome"] = {"type": "unsupported", "assembly": "hg38", "columns": {}}
    project.set_datasource_metadata(metadata)
    with pytest.raises(ValueError, match="Unsupported genome metadata type"):
        add_features_to_project(project, "ds", str(tmp_path / "missing.bed"))

    metadata["genome"] = {
        "type": "interval",
        "assembly": "hg38",
        "columns": {"chr": "chr", "start": "start"},
    }
    project.set_datasource_metadata(metadata)
    with pytest.raises(ValueError, match="missing required columns"):
        add_features_to_project(project, "ds", str(tmp_path / "missing.bed"))
