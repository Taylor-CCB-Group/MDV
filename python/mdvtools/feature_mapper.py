"""Genomic feature mapping helpers built on pyranges.

Primary API:
- map_sv_genes(sv_df, genes_bed_file, window_bp=0)
- add_features_to_project(project, datasource, genes_bed_file, ...)

Input dataframe requirements:
- index: used as sv_id
- columns: chr1, pos1, chr2, pos2, svtype

Returned dataframe:
- map_sv_genes(): columns sv_id, genes
"""

from __future__ import annotations

import re
from collections import defaultdict
from typing import Any

# Regex matching gene identifiers that are accession numbers rather than names,
# e.g. ENSG00000123456, ENST..., ENSP..., LRG_..., or plain numeric IDs.
_UNNAMED_GENE_RE = re.compile(
    r"^(ENS[A-Z]+\d+|LRG_\d+|\d+)$",
    re.IGNORECASE,
)

import pandas as pd

try:
    import pyranges as pr
except ModuleNotFoundError as exc:  # pragma: no cover
    raise ModuleNotFoundError(
        "pyranges is required for mdvtools.feature_mapper. "
        "Install with: pip install pyranges"
    ) from exc


def _chr_norm(chrom: Any) -> str:
    if pd.isna(chrom):
        return ""
    text = str(chrom).strip()
    if text == "":
        return ""
    core = text[3:] if text.lower().startswith("chr") else text
    return f"chr{core}"


def _to_int(value: Any) -> int | None:
    if pd.isna(value):
        return None
    text = str(value).strip()
    if text == "":
        return None
    try:
        return int(float(text))
    except ValueError:
        return None


def _read_genes_bed_df(genes_bed_file: str) -> pd.DataFrame:
    genes_df = pd.read_csv(
        genes_bed_file,
        sep="\t",
        header=None,
        comment="#",
        usecols=range(4),
        names=["Chromosome", "Start", "End", "gene_name"],
    )
    if genes_df.empty:
        return genes_df

    genes_df["Chromosome"] = genes_df["Chromosome"].map(_chr_norm)
    genes_df["Start"] = pd.to_numeric(genes_df["Start"], errors="coerce")
    genes_df["End"] = pd.to_numeric(genes_df["End"], errors="coerce")
    genes_df["gene_name"] = genes_df["gene_name"].astype(str).str.strip()

    valid = (
        genes_df["Chromosome"].ne("")
        & genes_df["Start"].notna()
        & genes_df["End"].notna()
        & genes_df["gene_name"].ne("")
    )
    genes_df = genes_df.loc[valid].copy()
    if genes_df.empty:
        return genes_df

    genes_df["Start"] = genes_df["Start"].astype(int)
    genes_df["End"] = genes_df["End"].astype(int)
    genes_df = genes_df.loc[genes_df["End"] > genes_df["Start"]].copy()
    return genes_df


def _build_interval_bed_df(interval_df: pd.DataFrame) -> pd.DataFrame:
    rows: list[list[Any]] = []
    required = ["chr", "start", "end"]
    missing = [c for c in required if c not in interval_df.columns]
    if missing:
        raise ValueError(f"Missing required interval dataframe columns: {missing}")

    for row_id, rec in interval_df.iterrows():
        chrom = _chr_norm(rec.get("chr"))
        start = _to_int(rec.get("start"))
        end = _to_int(rec.get("end"))
        if not chrom or start is None or end is None:
            continue
        bed_start = max(0, min(start, end))
        bed_end = max(start, end)
        if bed_end > bed_start:
            rows.append([chrom, bed_start, bed_end, str(row_id)])

    return pd.DataFrame(
        rows,
        columns=pd.Index(["chrom", "start", "end", "row_id"]),
    )


def _build_sv_bed_df(
    sv_df: pd.DataFrame, window_bp: int, inversion_as_bnd: bool = False
) -> pd.DataFrame:
    rows: list[list[Any]] = []
    required = ["chr1", "pos1", "chr2", "pos2", "svtype"]
    missing = [c for c in required if c not in sv_df.columns]
    if missing:
        raise ValueError(f"Missing required SV dataframe columns: {missing}")

    if window_bp < 0:
        raise ValueError("window_bp must be >= 0")

    for sv_id, rec in sv_df.iterrows():
        sv_id_str = str(sv_id)
        svtype = str(rec.get("svtype", "")).strip().upper()
        c1 = _chr_norm(rec.get("chr1"))
        c2 = _chr_norm(rec.get("chr2"))
        p1 = _to_int(rec.get("pos1"))
        p2 = _to_int(rec.get("pos2"))

        if svtype in {"BND", "TRA"} or (inversion_as_bnd and svtype == "INV"):
            if c1 and p1 is not None:
                start = max(0, p1 - 1 - window_bp)
                end = p1 + window_bp
                if end > start:
                    rows.append([c1, start, end, sv_id_str])
            if c2 and p2 is not None:
                start = max(0, p2 - 1 - window_bp)
                end = p2 + window_bp
                if end > start:
                    rows.append([c2, start, end, sv_id_str])
        else:
            if not c1 or p1 is None:
                continue
            right = p2 if p2 is not None else p1
            start = max(0, min(p1, right) - 1)
            end = max(p1, right)
            if end > start:
                rows.append([c1, start, end, sv_id_str])

    return pd.DataFrame(
        rows,
        columns=pd.Index(["chrom", "start", "end", "sv_id"]),
    )


def _map_rows_to_features(
    row_bed_df: pd.DataFrame,
    row_ids: pd.Index,
    genes_bed_file: str,
    row_id_column: str,
    output_column: str,
    filter_unnamed: bool = False,
) -> pd.DataFrame:
    normalized_ids = row_ids.map(str)
    if row_bed_df.empty:
        return pd.DataFrame(
            {row_id_column: normalized_ids, output_column: [""] * len(normalized_ids)}
        )

    genes_df = _read_genes_bed_df(genes_bed_file)
    if genes_df.empty:
        return pd.DataFrame(
            {row_id_column: normalized_ids, output_column: [""] * len(normalized_ids)}
        )

    row_ranges = pr.PyRanges(
        row_bed_df.rename(columns={"chrom": "Chromosome", "start": "Start", "end": "End"})
    )
    gene_ranges = pr.PyRanges(genes_df)
    intersections = row_ranges.join(gene_ranges).df

    features_by_row: defaultdict[str, set[str]] = defaultdict(set)
    for _, interval in intersections.iterrows():
        row_id = str(interval[row_id_column]).strip()
        gene_name = str(interval["gene_name"]).strip()
        if row_id and gene_name:
            if filter_unnamed and _UNNAMED_GENE_RE.match(gene_name):
                continue
            features_by_row[row_id].add(gene_name)

    out = pd.DataFrame({row_id_column: normalized_ids})
    out[output_column] = out[row_id_column].map(
        lambda rid: ",".join(sorted(features_by_row.get(rid, set())))
    )
    return out


def _get_required_genome_columns(genome_info: Any) -> tuple[str, list[str]]:
    if not isinstance(genome_info, dict):
        raise ValueError("Datasource genome metadata must be a dictionary")

    genome_type = genome_info.get("type")
    columns = genome_info.get("columns")
    if not isinstance(columns, dict):
        raise ValueError("Datasource genome metadata must include a columns mapping")

    if genome_type == "interval":
        required = ["chr", "start", "end"]
    elif genome_type == "sv":
        required = ["chr1", "pos1", "chr2", "pos2", "svtype", "length"]
    else:
        raise ValueError(f"Unsupported genome metadata type: {genome_type!r}")

    missing = [name for name in required if name not in columns]
    if missing:
        raise ValueError(f"Datasource genome metadata is missing required columns: {missing}")

    resolved_columns = [str(columns[name]) for name in required]
    return str(genome_type), resolved_columns


def add_features_to_project(
    project: Any,
    datasource: str,
    genes_bed_file: str,
    window_bp: int = 0,
    filter_unnamed: bool = True,
    inversion_as_bnd: bool = True,
    column_name: str = "features",
) -> None:
    """Add mapped genomic features to an MDVProject datasource.

    The datasource must provide genome metadata using the current schema in
    datasource["genome"], with type "interval" or "sv" and a "columns" mapping.
    """
    ds_metadata = project.get_datasource_metadata(datasource)
    genome_info = ds_metadata.get("genome")
    if genome_info is None:
        raise ValueError(f"Datasource {datasource!r} is missing genome metadata")

    genome_type, required_columns = _get_required_genome_columns(genome_info)
    source_df = project.get_datasource_as_dataframe(datasource, columns=required_columns)

    if genome_type == "interval":
        rename_map = dict(zip(required_columns, ["chr", "start", "end"]))
        row_bed_df = _build_interval_bed_df(source_df.rename(columns=rename_map))
        features_df = _map_rows_to_features(
            row_bed_df,
            source_df.index,
            genes_bed_file,
            row_id_column="row_id",
            output_column=column_name,
            filter_unnamed=filter_unnamed,
        )
    else:
        rename_map = dict(
            zip(required_columns, ["chr1", "pos1", "chr2", "pos2", "svtype", "length"])
        )
        sv_df = source_df.rename(columns=rename_map)
        row_bed_df = _build_sv_bed_df(
            sv_df, window_bp, inversion_as_bnd=inversion_as_bnd
        )
        features_df = _map_rows_to_features(
            row_bed_df,
            sv_df.index,
            genes_bed_file,
            row_id_column="sv_id",
            output_column=column_name,
            filter_unnamed=filter_unnamed,
        )

    project.set_column(
        datasource,
        {
            "name": column_name,
            "field": column_name,
            "datatype": "multitext",
            "delimiter": ",",
        },
        features_df[column_name].tolist(),
    )


def map_sv_genes(
    sv_df: pd.DataFrame,
    genes_bed_file: str,
    window_bp: int = 0,
    filter_unnamed: bool = True,
    inversion_as_bnd: bool = True,
) -> pd.DataFrame:
    """Map SVs to overlapping genes using pyranges.

    Parameters
    ----------
    sv_df:
        Input dataframe with index as sv_id and columns chr1,pos1,chr2,pos2,svtype.
    genes_bed_file:
        Path to 4-column BED file (chr, start, end, gene_name).
    window_bp:
        Flank size around BND/TRA breakpoints.
    filter_unnamed:
        When True, discard gene names that look like accession IDs rather than
        human-readable names (e.g. ENSG00000123456, ENST..., LRG_...).
    inversion_as_bnd:
        When True, treat INV like BND/TRA (emit two breakpoints, not a range; SV is in gene if either end overlaps).

    Returns
    -------
    pandas.DataFrame
        Two columns: sv_id and genes (comma-delimited unique genes).
    """
    sv_bed_df = _build_sv_bed_df(sv_df, window_bp, inversion_as_bnd=inversion_as_bnd)
    return _map_rows_to_features(
        sv_bed_df,
        sv_df.index,
        genes_bed_file,
        row_id_column="sv_id",
        output_column="genes",
        filter_unnamed=filter_unnamed,
    )
