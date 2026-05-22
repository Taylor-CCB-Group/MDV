"""SV to gene mapping helpers built on pybedtools.

Primary API:
- map_sv_genes(sv_df, genes_bed_file, window_bp=0)

Input dataframe requirements:
- index: used as sv_id
- columns: chr1, pos1, chr2, pos2, svtype

Returned dataframe:
- columns: sv_id, genes
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
    import pybedtools
except ModuleNotFoundError as exc:  # pragma: no cover
    raise ModuleNotFoundError(
        "pybedtools is required for scripts.sv_gene_mapper. "
        "Install with: pip install pybedtools"
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


def _build_sv_bed_df(sv_df: pd.DataFrame, window_bp: int, inversion_as_bnd: bool = False) -> pd.DataFrame:
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


def map_sv_genes(
    sv_df: pd.DataFrame,
    genes_bed_file: str,
    window_bp: int = 0,
    filter_unnamed: bool = False,
    inversion_as_bnd: bool = False,
) -> pd.DataFrame:
    """Map SVs to overlapping genes using pybedtools.

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
    all_sv_ids = sv_df.index.map(str)

    if sv_bed_df.empty:
        return pd.DataFrame({"sv_id": all_sv_ids, "genes": [""] * len(all_sv_ids)})

    sv_bt = pybedtools.BedTool.from_dataframe(sv_bed_df)
    genes_bt = pybedtools.BedTool(genes_bed_file)
    intersections = sv_bt.intersect(genes_bt, wa=True, wb=True)

    genes_by_sv: defaultdict[str, set[str]] = defaultdict(set)
    for interval in intersections:
        fields = interval.fields
        if len(fields) < 8:
            continue
        sv_id = fields[3].strip()
        gene_name = fields[7].strip()
        if sv_id and gene_name:
            if filter_unnamed and _UNNAMED_GENE_RE.match(gene_name):
                continue
            genes_by_sv[sv_id].add(gene_name)

    out = pd.DataFrame({"sv_id": all_sv_ids})
    out["genes"] = out["sv_id"].map(
        lambda sid: ",".join(sorted(genes_by_sv.get(sid, set())))
    )
    return out

