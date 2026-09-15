"""Infer experimental-design fields from an MDV project and add template views.

Usage:
  mdvtools create-default-views /path/to/mdv_project
  python -m mdvtools.template_views --project /path/to/mdv_project

See python/mdvtools/template_views/README.md for what is written and how.
"""

from __future__ import annotations

import argparse
import json
import random
import shutil
import string
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

VIEW_NAMES = [
    "1 Study overview",
    "2 Cell atlas",
    "3 Marker genes and signatures",
    "4 Tissue and disease context",
    "5 RNA-protein concordance",
]

MARKER_PANEL = [
    "CD14",
    "FCGR3A",
    "CD3D",
    "CD3E",
    "MS4A1",
    "CD19",
    "NKG7",
    "GNLY",
    "CD1C",
    "CLEC9A",
    "C1QA",
    "C1QB",
    "LYZ",
    "S100A8",
    "S100A9",
    "VCAN",
    "FCN1",
    "MERTK",
    "LYVE1",
    "APOE",
    "TREM2",
    "IL3RA",
    "KIT",
    "TNF",
    "IL1B",
    "ISG15",
    "CD4",
    "CD8A",
    "FOXP3",
    "MKI67",
]
FEATURE_UMAP_GENES = ["CD14", "FCGR3A", "CD3D", "MS4A1", "CD1C"]
PROTEIN_PANEL = [
    "adt_CD14",
    "CD14",
    "adt_CD16",
    "adt_CD3",
    "adt_CD19",
    "adt_CD4",
    "adt_CD8",
    "adt_CD56",
    "adt_HLA-DR",
    "adt_CD1c",
    "adt_CD11c",
    "adt_CD86",
]

CELL_TYPE_PATTERNS = [
    "annotation_caf",
    "annotation",
    "cell_type",
    "celltype",
    "leiden_res1",
    "leiden",
    "cluster",
    "final_analysis",
]
BROAD_TYPE_PATTERNS = ["sub_bucket_caf", "sub_bucket", "major", "leiden_res1"]
TISSUE_PATTERNS = ["tissue_simple", "tissue", "Run_Tissue_name", "spatial_region"]
DISEASE_PATTERNS = ["diagnosis", "disease"]
TREATMENT_PATTERNS = ["treatment_simple", "treatment"]
RESPONSE_PATTERNS = ["response"]
INFLAMMATION_PATTERNS = ["inflammation_status", "inflammation"]
WORKSTREAM_PATTERNS = ["workstream"]
SEX_PATTERNS = ["gender", "sex"]
SAMPLE_PATTERNS = ["sample_id", "cart_id", "Run_Tissue_name", "slide_ID"]
QC_PATTERNS = [
    "n_genes_by_counts",
    "pct_counts_mt",
    "total_counts",
    "nCount_RNA",
    "nFeature_RNA",
    "nCount_negprobes",
]
FACTOR_PATTERNS = [
    "treatment_simple",
    "treatment",
    "disease_grp_treatment",
    "diagnosis",
    "disease",
    "inflammation_status",
    "inflammation",
    "response",
    "tissue_simple",
    "tissue",
    "Run_Tissue_name",
    "spatial_region",
]
CATEGORICAL_DTYPES = {"text", "text16", "multitext"}
NUMERIC_DTYPES = {"integer", "double", "int32"}
MARKER_DS_NAME = "cluster_markers"
MARKER_CACHE = "cluster_markers.json"
MARKERS_PER_CLUSTER = 20
VARYING_GENES_PER_FACTOR = 15
FEATURE_PLOT_CAP = 24
MISSING_LABELS = {"nd", "nan", "na", "", "undetermined or na", "null"}
SKIP_GENE_PREFIXES = ("MT-", "RPL", "RPS")
SKIP_GENES = {"MALAT1", "NEAT1"}
SEX_GENES = {
    "XIST",
    "DDX3Y",
    "RPS4Y1",
    "EIF1AY",
    "UTY",
    "KDM5D",
    "ZFY",
    "TMSB4Y",
    "NLGN4Y",
    "USP9Y",
}
MATRIX_PREFERENCE = ["rna_logged_counts", "rna_raw_counts", "gs", "rna_expr"]


def _id() -> str:
    return "".join(random.choices(string.ascii_letters, k=6))


def _chart(
    chart_type: str,
    title: str,
    param: Any,
    gsposition: list[int],
    gssize: list[int],
    **extra: Any,
) -> dict[str, Any]:
    data: dict[str, Any] = {
        "type": chart_type,
        "title": title,
        "param": param,
        "size": [400, 300],
        "position": [0, 0],
        "gsposition": gsposition,
        "gssize": gssize,
        "id": _id(),
    }
    data.update(extra)
    return data


def axis() -> dict[str, Any]:
    return {"label": "", "textSize": 13, "tickfont": 10}


def scatter_axis(label: str, size: int = 30) -> dict[str, Any]:
    return {
        "size": size,
        "label": label,
        "textSize": 13,
        "textsize": 13,
        "tickfont": 10,
    }


def row(title: str, field: str, pos: list[int], size: list[int]) -> dict[str, Any]:
    return _chart("row_chart", title, [field], pos, size, axis={"x": axis()})


def stacked(title: str, params: list[str], pos: list[int], size: list[int]) -> dict[str, Any]:
    return _chart(
        "stacked_row_chart",
        title,
        params,
        pos,
        size,
        axis={"x": axis(), "y": axis()},
        color_legend={"display": True, "pos": [8, 8]},
    )


def umap(
    title: str,
    color_by: str,
    pos: list[int],
    size: list[int],
    coords: list[str],
) -> dict[str, Any]:
    return _chart(
        "wgl_scatter_plot",
        title,
        coords,
        pos,
        size,
        color_by=color_by,
        default_color="#377eb8",
        axis={"x": scatter_axis(coords[0]), "y": scatter_axis(coords[1], 45)},
        color_legend={"display": True, "pos": [12, 12]},
        radius=3,
        opacity=0.8,
        brush="poly",
    )


def dot(title: str, params: list[Any], pos: list[int], size: list[int]) -> dict[str, Any]:
    return _chart(
        "dot_plot",
        title,
        params,
        pos,
        size,
        axis={"x": axis(), "y": axis()},
        color_scale={"log": False},
        color_legend={"display": True, "pos": [40, 10]},
        fraction_legend={"display": True, "pos": [0, 0]},
    )


def heatmap(title: str, params: list[Any], pos: list[int], size: list[int]) -> dict[str, Any]:
    return _chart("heat_map", title, params, pos, size)


def violin(title: str, params: list[Any], pos: list[int], size: list[int]) -> dict[str, Any]:
    return _chart(
        "violin_plot",
        title,
        params,
        pos,
        size,
        box_visible=True,
        meanline_visible=True,
    )


def histogram(
    title: str,
    field: str,
    display_min: float,
    display_max: float,
    pos: list[int],
    size: list[int],
) -> dict[str, Any]:
    return _chart(
        "bar_chart",
        title,
        [field],
        pos,
        size,
        bin_number=50,
        display_min=display_min,
        display_max=display_max,
        x_axis={"size": 30, "label": field, "textSize": 13, "tickfont": 10},
        y_axis={
            "size": 45,
            "label": "cells",
            "textSize": 13,
            "tickfont": 10,
            "rotate_labels": False,
        },
    )


def selection(title: str, params: list[str], pos: list[int], size: list[int]) -> dict[str, Any]:
    return _chart("selection_dialog", title, params, pos, size)


def table(title: str, params: list[str], pos: list[int], size: list[int]) -> dict[str, Any]:
    return _chart("table_chart", title, params, pos, size, include_index=False)


def textbox(title: str, text: str, pos: list[int], size: list[int]) -> dict[str, Any]:
    return _chart("text_box_chart", title, [], pos, size, text=text)


def abundance(title: str, params: list[str], pos: list[int], size: list[int]) -> dict[str, Any]:
    return _chart(
        "custom_box_plot",
        title,
        params,
        pos,
        size,
        x_axis={"labels": [""], "title": params[0] if params else ""},
        y_axis={"labels": ["Abundance"], "title": "fraction"},
    )


def load_json(path: Path) -> Any:
    return json.loads(path.read_text())


def save_json(path: Path, data: Any) -> None:
    path.write_text(json.dumps(data, indent=2) + "\n")


def suffix(field: str) -> str:
    return field.split(":")[-1]


def prefix(field: str) -> str | None:
    if ":" in field:
        return field.split(":", 1)[0]
    return None


@dataclass
class LinkInfo:
    ds_name: str
    name_column: str
    subgroup: str
    names: list[str]
    matrix: str | None = None
    matrix_h5: str | None = None


@dataclass
class MarkerResult:
    rows: list[dict[str, Any]] = field(default_factory=list)
    factor_genes: dict[str, list[str]] = field(default_factory=dict)
    table_fields: list[str] = field(default_factory=list)


def unique_top_markers(
    markers: MarkerResult | None,
    per_cluster: int = 2,
    cap: int = FEATURE_PLOT_CAP,
) -> list[tuple[str, str]]:
    """Unique rank-1..N cluster-vs-rest genes, first cluster wins, capped."""
    if not markers:
        return []
    by_cluster: dict[str, list[tuple[int, str]]] = {}
    for row in markers.rows:
        if row.get("contrast") != "cluster vs rest":
            continue
        gene = str(row.get("gene") or "")
        cluster = str(row.get("cluster") or "")
        if not gene or not cluster:
            continue
        try:
            rank = int(row.get("rank") or 0)
        except (TypeError, ValueError):
            continue
        if rank < 1 or rank > per_cluster:
            continue
        by_cluster.setdefault(cluster, []).append((rank, gene))
    picked: list[tuple[str, str]] = []
    seen: set[str] = set()
    for cluster, hits in by_cluster.items():
        for _rank, gene in sorted(hits, key=lambda t: t[0]):
            if gene in seen:
                continue
            seen.add(gene)
            picked.append((gene, cluster))
            if len(picked) >= cap:
                return picked
    return picked


@dataclass
class Roles:
    obs: str
    columns: dict[str, dict[str, Any]]
    fields: list[str]
    rna: LinkInfo | None = None
    protein: LinkInfo | None = None
    cell_type: str | None = None
    broad_type: str | None = None
    tissue: str | None = None
    disease: str | None = None
    treatment: str | None = None
    response: str | None = None
    inflammation: str | None = None
    workstream: str | None = None
    sex: str | None = None
    sample_id: str | None = None
    embedding: list[str] | None = None
    protein_embedding: list[str] | None = None
    qc: list[str] = field(default_factory=list)
    scores: list[str] = field(default_factory=list)
    factors: list[str] = field(default_factory=list)
    gene_wrappers: dict[str, str] = field(default_factory=dict)
    protein_wrappers: dict[str, str] = field(default_factory=dict)
    markers: MarkerResult | None = None
    skipped: list[str] = field(default_factory=list)
    n_cells: int | None = None


def ds_by_name(datasources: list[dict[str, Any]], name: str) -> dict[str, Any] | None:
    for d in datasources:
        if d.get("name") == name:
            return d
    return None


def name_values(ds: dict[str, Any], name_column: str) -> list[str]:
    for col in ds.get("columns") or []:
        if col.get("field") == name_column:
            values = col.get("values")
            if isinstance(values, list):
                return [str(v) for v in values]
    return []


def n_values(col: dict[str, Any]) -> int | None:
    values = col.get("values")
    if isinstance(values, list):
        return len(values)
    return None


def minmax(col: dict[str, Any]) -> tuple[float, float] | None:
    mm = col.get("minMax")
    if isinstance(mm, list) and len(mm) == 2:
        return float(mm[0]), float(mm[1])
    return None


def pick_obs(datasources: list[dict[str, Any]]) -> dict[str, Any]:
    named = ds_by_name(datasources, "cells")
    if named:
        return named
    for d in datasources:
        links = d.get("links") or {}
        if any(
            isinstance(v, dict) and "rows_as_columns" in v for v in links.values()
        ):
            return d
    if not datasources:
        raise ValueError("Project has no datasources")
    return datasources[0]


def pick_subgroup(subgroups: dict[str, Any], preferred: list[str]) -> str | None:
    if not subgroups:
        return None
    for key in preferred:
        if key in subgroups:
            return key
    return next(iter(subgroups))


def match_field(
    fields: list[str],
    patterns: list[str],
    prefer_prefix: str = "rna",
) -> str | None:
    for pat in patterns:
        exact = [f for f in fields if f == pat]
        suff = [f for f in fields if suffix(f) == pat]
        hits = exact + [f for f in suff if f not in exact]
        if not hits:
            continue
        pref = [f for f in hits if prefix(f) == prefer_prefix]
        return (pref or hits)[0]
    return None


def pick_nbclust_field(fields: list[str], cols: dict[str, dict[str, Any]]) -> str | None:
    cands = [f for f in fields if "RNA_nbclust" in f and f.endswith("_clusters")]
    if not cands:
        return None

    def score(fid: str) -> tuple[int, int]:
        vals = [
            str(v)
            for v in (cols.get(fid) or {}).get("values") or []
            if not is_missing_label(str(v))
        ]
        named = sum(1 for v in vals if len(v) > 1)
        return (named, len(vals))

    return max(cands, key=score)


def pick_spatial_niche(fields: list[str]) -> str | None:
    cands = [f for f in fields if f.startswith("spatialclust_") and f.endswith("_assignments")]
    return cands[0] if cands else None


def match_embedding_pair(
    fields: list[str],
    prefer_prefix: str | None,
    preferred_suffixes: list[tuple[str, str]],
) -> list[str] | None:
    by_suffix = {suffix(f): f for f in fields}

    def candidates(s1: str, s2: str) -> list[list[str]]:
        out: list[list[str]] = []
        ones = [f for f in fields if suffix(f) == s1]
        for a in ones:
            b_suffix = s2
            pref = prefix(a)
            expected = f"{pref}:{b_suffix}" if pref else b_suffix
            if expected in fields:
                out.append([a, expected])
            elif b_suffix in by_suffix and prefix(by_suffix[b_suffix]) == pref:
                out.append([a, by_suffix[b_suffix]])
        return out

    for s1, s2 in preferred_suffixes:
        pairs = candidates(s1, s2)
        if prefer_prefix:
            pref = [p for p in pairs if prefix(p[0]) == prefer_prefix]
            if pref:
                return pref[0]
        if pairs:
            return pairs[0]

    # any umap *_1 / *_2
    ones = [f for f in fields if "umap" in suffix(f).lower() and suffix(f).endswith("_1")]
    for a in ones:
        s2 = suffix(a)[:-1] + "2"
        expected = f"{prefix(a)}:{s2}" if prefix(a) else s2
        if expected in fields:
            if prefer_prefix and prefix(a) != prefer_prefix:
                continue
            return [a, expected]
    if ones:
        a = ones[0]
        s2 = suffix(a)[:-1] + "2"
        expected = f"{prefix(a)}:{s2}" if prefix(a) else s2
        if expected in fields:
            return [a, expected]
    return None


def wrapper(subgroup: str, name: str, index: int) -> str:
    return f"{subgroup}|{name}({subgroup})|{index}"


def resolve_named_wrappers(link: LinkInfo, genes: list[str]) -> dict[str, str]:
    index = {n: i for i, n in enumerate(link.names)}
    found: dict[str, str] = {}
    for gene in genes:
        if gene in index and gene not in found:
            found[gene] = wrapper(link.subgroup, gene, index[gene])
    return found


def resolve_gene_wrappers(link: LinkInfo) -> dict[str, str]:
    return resolve_named_wrappers(link, MARKER_PANEL)


def resolve_protein_wrappers(link: LinkInfo) -> dict[str, str]:
    index = {n: i for i, n in enumerate(link.names)}
    lower = {n.lower(): n for n in link.names}
    found: dict[str, str] = {}
    for token in PROTEIN_PANEL:
        name = None
        if token in index:
            name = token
        else:
            adt = token if token.startswith("adt_") else f"adt_{token}"
            bare = token[4:] if token.startswith("adt_") else token
            for cand in (token, adt, bare, adt.lower(), bare.lower()):
                if cand in index:
                    name = cand
                    break
                if cand.lower() in lower:
                    name = lower[cand.lower()]
                    break
        if name and name not in found.values():
            found[name] = wrapper(link.subgroup, name, index[name])
    if not found and link.names:
        name = link.names[0]
        found[name] = wrapper(link.subgroup, name, 0)
    return found


def is_missing_label(label: str) -> bool:
    return str(label).strip().lower() in MISSING_LABELS


def skip_housekeeping(name: str) -> bool:
    upper = name.upper()
    if upper in SKIP_GENES:
        return True
    if name.startswith(("SystemControl", "Negative", "FalseCode", "NegPrb")):
        return True
    return any(upper.startswith(p) for p in SKIP_GENE_PREFIXES)


def skip_factor_gene(name: str) -> bool:
    if name.upper() in SEX_GENES:
        return True
    return name.startswith("FP") and len(name) > 2 and name[2].isdigit()


def read_codes(h5: Any, obs: str, field: str) -> Any:
    import numpy as np

    if field in h5[obs]:
        return np.asarray(h5[obs][field], dtype=np.int32)
    raise KeyError(f"{obs}/{field} not in datafile.h5")


def encode_groups(values: list[str], codes: Any) -> tuple[Any, list[str]]:
    import numpy as np

    usable: list[str] = []
    remap = np.full(max(len(values), 1), -1, dtype=np.int32)
    for i, lab in enumerate(values):
        if is_missing_label(str(lab)):
            continue
        remap[i] = len(usable)
        usable.append(str(lab))
    raw = np.asarray(codes, dtype=np.int32)
    out = np.full(raw.shape[0], -1, dtype=np.int32)
    valid = (raw >= 0) & (raw < len(remap))
    out[valid] = remap[raw[valid]]
    return out, usable


def iter_sparse_columns(grp: Any, chunk_nnz: int = 8_000_000):
    import numpy as np

    p = np.asarray(grp["p"], dtype=np.int64)
    n_genes = len(p) - 1
    gene = 0
    while gene < n_genes:
        start = int(p[gene])
        end_gene = gene + 1
        while end_gene < n_genes and int(p[end_gene + 1]) - start <= chunk_nnz:
            end_gene += 1
        stop = int(p[end_gene])
        i = np.asarray(grp["i"][start:stop], dtype=np.int64)
        x = np.asarray(grp["x"][start:stop], dtype=np.float32)
        for g in range(gene, end_gene):
            a = int(p[g]) - start
            b = int(p[g + 1]) - start
            yield g, i[a:b], x[a:b]
        gene = end_gene


def accumulate_groups(codes: Any, n_groups: int, rows: Any, vals: Any) -> tuple[Any, Any]:
    import numpy as np

    group = codes[rows]
    keep = group >= 0
    if not np.any(keep):
        zeros = np.zeros(n_groups, dtype=np.float64)
        return zeros, np.zeros(n_groups, dtype=np.int64)
    group = group[keep]
    vals = vals[keep]
    sums = np.bincount(group, weights=vals.astype(np.float64), minlength=n_groups)
    nnz = np.bincount(group, minlength=n_groups)
    return sums, nnz


def compute_marker_tables(
    project_path: Path,
    roles: Roles,
) -> MarkerResult:
    import heapq

    import h5py
    import numpy as np

    result = MarkerResult()
    if not roles.rna or not roles.rna.matrix or not roles.cell_type:
        roles.skipped.append("cluster markers (no RNA matrix or cell type)")
        return result
    col = roles.columns.get(roles.cell_type) or {}
    cluster_values = col.get("values")
    if not isinstance(cluster_values, list) or len(cluster_values) < 2:
        roles.skipped.append("cluster markers (cell type has no values)")
        return result

    h5_path = project_path / "datafile.h5"
    if not h5_path.is_file():
        roles.skipped.append("cluster markers (no datafile.h5)")
        return result

    names = roles.rna.names
    matrix_h5 = roles.rna.matrix_h5 or roles.rna.matrix
    print(f"Computing markers from {roles.obs}/{matrix_h5} ({len(names)} genes)...")
    rows: list[dict[str, Any]] = []
    factor_genes: dict[str, list[str]] = {}

    with h5py.File(h5_path, "r") as h5:
        if roles.obs not in h5:
            roles.skipped.append(f"cluster markers (missing {roles.obs})")
            return result
        obs_obj = h5[roles.obs]
        if not isinstance(obs_obj, h5py.Group):
            roles.skipped.append(f"cluster markers ({roles.obs} is not a group)")
            return result
        obs_grp = obs_obj
        if matrix_h5 not in obs_grp:
            fallback = next(
                (c for c in (roles.rna.matrix, "gene_scores", "rna_logged_counts", "gs") if c and c in obs_grp),
                None,
            )
            if not fallback:
                roles.skipped.append(f"cluster markers (missing {matrix_h5})")
                return result
            matrix_h5 = fallback
        if roles.rna.matrix:
            roles.rna.matrix_h5 = matrix_h5
        cluster_codes_raw = read_codes(h5, roles.obs, roles.cell_type)
        cluster_codes, clusters = encode_groups([str(v) for v in cluster_values], cluster_codes_raw)
        n_clusters = len(clusters)
        if n_clusters < 2:
            roles.skipped.append("cluster markers (fewer than 2 usable clusters)")
            return result
        cluster_n = np.bincount(cluster_codes[cluster_codes >= 0], minlength=n_clusters).astype(np.int64)
        n_cells = int(cluster_codes.shape[0])

        factor_specs: list[tuple[str, Any, list[str], Any]] = []
        for fid in roles.factors:
            fcol = roles.columns.get(fid) or {}
            fvals = fcol.get("values")
            if not isinstance(fvals, list):
                continue
            try:
                raw = read_codes(h5, roles.obs, fid)
            except KeyError:
                continue
            codes, labels = encode_groups([str(v) for v in fvals], raw)
            if len(labels) < 2:
                continue
            combo = np.full(n_cells, -1, dtype=np.int32)
            ok = (cluster_codes >= 0) & (codes >= 0)
            combo[ok] = cluster_codes[ok] * len(labels) + codes[ok]
            combo_n = np.bincount(combo[combo >= 0], minlength=n_clusters * len(labels)).astype(np.int64)
            factor_specs.append((fid, combo, labels, combo_n))
            print(f"  factor {fid}: {labels}")

        print(f"  clusters ({n_clusters}): {clusters}")
        cluster_best: list[list[tuple[float, int, tuple[float, float, float, float]]]] = [
            [] for _ in clusters
        ]
        factor_best: dict[str, list[tuple[float, int]]] = {fid: [] for fid, *_ in factor_specs}
        matrix_obj = obs_grp[matrix_h5]
        if not isinstance(matrix_obj, h5py.Group):
            roles.skipped.append(f"cluster markers ({matrix_h5} is not a group)")
            return result
        grp = matrix_obj
        n_genes = len(np.asarray(grp["p"])) - 1
        reported = 0
        for gene_i, cell_i, vals in iter_sparse_columns(grp):
            if gene_i >= len(names):
                break
            gene = names[gene_i]
            if skip_housekeeping(gene) or cell_i.size < 50:
                continue
            sums, nnz = accumulate_groups(cluster_codes, n_clusters, cell_i, vals)
            total_sum = float(sums.sum())
            total_nnz = int(nnz.sum())
            for k, cluster in enumerate(clusters):
                n_in = int(cluster_n[k])
                n_out = n_cells - n_in
                if n_in < 20 or n_out < 20:
                    continue
                mean_in = float(sums[k] / max(n_in, 1))
                mean_out = float((total_sum - sums[k]) / max(n_out, 1))
                pct_in = float(nnz[k] / max(n_in, 1))
                pct_out = float((total_nnz - nnz[k]) / max(n_out, 1))
                if pct_in < 0.1 or mean_in <= mean_out:
                    continue
                score = (mean_in - mean_out) * (pct_in - pct_out + 0.05)
                item = (score, gene_i, (mean_in, mean_out, pct_in, pct_out))
                heap = cluster_best[k]
                if len(heap) < MARKERS_PER_CLUSTER:
                    heapq.heappush(heap, item)
                elif score > heap[0][0]:
                    heapq.heapreplace(heap, item)
            if not skip_factor_gene(gene):
                for fid, combo, labels, combo_n in factor_specs:
                    fsums, _fnz = accumulate_groups(combo, n_clusters * len(labels), cell_i, vals)
                    ranges: list[float] = []
                    n_lab = len(labels)
                    for k in range(n_clusters):
                        means: list[float] = []
                        for t in range(n_lab):
                            idx = k * n_lab + t
                            n = int(combo_n[idx])
                            if n < 40:
                                continue
                            means.append(float(fsums[idx] / n))
                        if len(means) >= 2:
                            ranges.append(max(means) - min(means))
                    if ranges:
                        fscore = float(sum(ranges) / len(ranges))
                        fheap = factor_best[fid]
                        fitem = (fscore, gene_i)
                        if len(fheap) < VARYING_GENES_PER_FACTOR:
                            heapq.heappush(fheap, fitem)
                        elif fscore > fheap[0][0]:
                            heapq.heapreplace(fheap, fitem)
            reported += 1
            if reported % 4000 == 0:
                print(f"  scored {gene_i + 1}/{n_genes} genes")

    for k, cluster in enumerate(clusters):
        top = sorted(cluster_best[k], key=lambda t: t[0], reverse=True)
        for rank, (score, gene_i, stats) in enumerate(top, start=1):
            rows.append(
                {
                    "contrast": "cluster vs rest",
                    "cluster": cluster,
                    "gene": names[gene_i],
                    "rank": rank,
                    "score": float(score),
                    "mean_in": stats[0],
                    "mean_out": stats[1],
                    "pct_in": stats[2],
                    "pct_out": stats[3],
                }
            )

    for fid, hits in factor_best.items():
        top = sorted(hits, key=lambda t: t[0], reverse=True)
        genes = [names[i] for _s, i in top]
        factor_genes[fid] = genes
        for rank, (score, gene_i) in enumerate(top, start=1):
            rows.append(
                {
                    "contrast": suffix(fid),
                    "cluster": "within cell type",
                    "gene": names[gene_i],
                    "rank": rank,
                    "score": float(score),
                    "mean_in": float("nan"),
                    "mean_out": float("nan"),
                    "pct_in": float("nan"),
                    "pct_out": float("nan"),
                }
            )

    result.rows = rows
    result.factor_genes = factor_genes
    result.table_fields = [
        "contrast",
        "cluster",
        "gene",
        "rank",
        "score",
        "mean_in",
        "mean_out",
        "pct_in",
        "pct_out",
    ]
    print(f"  {len(rows)} marker rows; factor genes: { {k: len(v) for k, v in factor_genes.items()} }")
    return result


def marker_result_to_json(result: MarkerResult) -> dict[str, Any]:
    return {
        "rows": result.rows,
        "factor_genes": result.factor_genes,
        "table_fields": result.table_fields,
    }


def marker_result_from_json(data: dict[str, Any]) -> MarkerResult:
    return MarkerResult(
        rows=list(data.get("rows") or []),
        factor_genes={str(k): list(v) for k, v in (data.get("factor_genes") or {}).items()},
        table_fields=list(data.get("table_fields") or []),
    )


def write_table_datasource(project_path: Path, name: str, rows: list[dict[str, Any]]) -> list[str]:
    import numpy as np
    import h5py

    if not rows:
        return []
    fields = list(rows[0].keys())
    columns: list[dict[str, Any]] = []
    h5_path = project_path / "datafile.h5"
    with h5py.File(h5_path, "a") as h5:
        if name in h5:
            del h5[name]
        gr = h5.create_group(name)
        for field in fields:
            vals = [r[field] for r in rows]
            if field == "rank":
                arr = np.asarray(vals, dtype=np.int32)
                gr.create_dataset(field, data=arr, dtype=np.int32)
                columns.append(
                    {
                        "datatype": "int32",
                        "name": field,
                        "field": field,
                        "minMax": [int(arr.min()), int(arr.max())],
                    }
                )
                continue
            if all(isinstance(v, (int, float)) and not isinstance(v, bool) for v in vals):
                arr = np.asarray(vals, dtype=np.float32)
                finite = arr[np.isfinite(arr)]
                mm = [float(finite.min()), float(finite.max())] if finite.size else [0.0, 0.0]
                gr.create_dataset(field, data=arr, dtype=np.float32)
                columns.append(
                    {
                        "datatype": "double",
                        "name": field,
                        "field": field,
                        "minMax": mm,
                        "quantiles": {},
                    }
                )
                continue
            strs = [str(v) for v in vals]
            uniq: list[str] = []
            seen: set[str] = set()
            for s in strs:
                if s not in seen:
                    seen.add(s)
                    uniq.append(s)
            vdict = {k: i for i, k in enumerate(uniq)}
            use16 = len(uniq) > 256
            codes = np.asarray([vdict[s] for s in strs], dtype=np.uint16 if use16 else np.uint8)
            gr.create_dataset(field, data=codes, dtype=codes.dtype)
            columns.append(
                {
                    "datatype": "text16" if use16 else "text",
                    "name": field,
                    "field": field,
                    "values": uniq,
                }
            )
    ds_path = project_path / "datasources.json"
    datasources = load_json(ds_path)
    if not isinstance(datasources, list):
        raise TypeError("datasources.json must be a list")
    entry = {"name": name, "columns": columns, "size": len(rows)}
    replaced = False
    for i, ds in enumerate(datasources):
        if ds.get("name") == name:
            datasources[i] = entry
            replaced = True
            break
    if not replaced:
        datasources.append(entry)
    save_json(ds_path, datasources)
    return [c["field"] for c in columns]


def load_or_compute_markers(
    project_path: Path,
    roles: Roles,
    skip: bool,
    recompute: bool,
) -> MarkerResult | None:
    cache = project_path / MARKER_CACHE
    if skip:
        return None
    if cache.is_file() and not recompute:
        print(f"Reusing marker cache {cache}")
        result = marker_result_from_json(load_json(cache))
        if result.rows:
            result.table_fields = write_table_datasource(project_path, MARKER_DS_NAME, result.rows) or result.table_fields
        return result
    result = compute_marker_tables(project_path, roles)
    if not result.rows:
        return result
    save_json(cache, marker_result_to_json(result))
    result.table_fields = write_table_datasource(project_path, MARKER_DS_NAME, result.rows)
    return result


def infer_roles(datasources: list[dict[str, Any]]) -> Roles:
    obs_ds = pick_obs(datasources)
    cols = {c.get("field"): c for c in obs_ds.get("columns") or [] if c.get("field")}
    fields = list(cols)
    roles = Roles(
        obs=str(obs_ds.get("name")),
        columns=cols,
        fields=fields,
        n_cells=obs_ds.get("size"),
    )

    links = obs_ds.get("links") or {}
    for target, spec in links.items():
        rac = (spec or {}).get("rows_as_columns") or {}
        subgroups = rac.get("subgroups") or {}
        if not subgroups:
            continue
        name_col = str(rac.get("name_column") or "name")
        target_ds = ds_by_name(datasources, str(target))
        names = name_values(target_ds, name_col) if target_ds else []
        tlow = str(target).lower()
        if roles.rna is None and any(k in tlow for k in ("rna", "gene")):
            sg = pick_subgroup(subgroups, ["rna_expr", "gs", "gene_scores"])
            matrix = pick_subgroup(subgroups, MATRIX_PREFERENCE)
            matrix_h5 = None
            if matrix:
                info = subgroups.get(matrix) or {}
                matrix_h5 = str(info.get("name") or matrix)
            if sg:
                roles.rna = LinkInfo(str(target), name_col, sg, names, matrix, matrix_h5)
        elif roles.protein is None and any(k in tlow for k in ("prot", "protein", "adt")):
            sg = pick_subgroup(subgroups, ["prot_expr", "prot_dsb", "pro_rc", "protein_raw_counts"])
            if sg:
                roles.protein = LinkInfo(str(target), name_col, sg, names)

    roles.cell_type = match_field(fields, CELL_TYPE_PATTERNS) or pick_nbclust_field(fields, cols)
    roles.broad_type = match_field(fields, BROAD_TYPE_PATTERNS) or pick_spatial_niche(fields)
    if roles.broad_type and roles.broad_type == roles.cell_type:
        roles.broad_type = match_field(
            [f for f in fields if f != roles.cell_type], BROAD_TYPE_PATTERNS
        )
    roles.tissue = match_field(fields, TISSUE_PATTERNS)
    roles.disease = match_field(fields, DISEASE_PATTERNS)
    roles.treatment = match_field(fields, TREATMENT_PATTERNS)
    roles.response = match_field(fields, RESPONSE_PATTERNS)
    roles.inflammation = match_field(fields, INFLAMMATION_PATTERNS)
    roles.workstream = match_field(fields, WORKSTREAM_PATTERNS)
    roles.sex = match_field(fields, SEX_PATTERNS)
    roles.sample_id = match_field(fields, SAMPLE_PATTERNS)
    seen_factors: set[str] = set()
    for pat in FACTOR_PATTERNS:
        fid = match_field(fields, [pat])
        if not fid or fid in seen_factors or fid == roles.cell_type:
            continue
        col = cols.get(fid) or {}
        if col.get("datatype") not in CATEGORICAL_DTYPES:
            continue
        nv = n_values(col)
        if nv is not None and (nv < 2 or nv > 20):
            continue
        roles.factors.append(fid)
        seen_factors.add(fid)

    roles.embedding = match_embedding_pair(
        fields,
        "rna",
        [
            ("X_umap_mindist_0.25_1", "X_umap_mindist_0.25_2"),
            ("X_umap_1", "X_umap_2"),
        ],
    )
    roles.protein_embedding = match_embedding_pair(
        fields,
        "prot",
        [("X_umap_1", "X_umap_2")],
    )
    if not roles.embedding:
        for a, b in (("x", "y"), ("spatial_1", "spatial_2"), ("global_1", "global_2")):
            if a in fields and b in fields:
                roles.embedding = [a, b]
                break

    for pat in QC_PATTERNS:
        fid = match_field(fields, [pat])
        if fid and minmax(cols[fid]):
            roles.qc.append(fid)

    for fid in fields:
        s = suffix(fid)
        if s.endswith("_score") and s != "doublet_scores":
            if cols[fid].get("datatype") in NUMERIC_DTYPES:
                if prefix(fid) in (None, "rna") or not any(
                    suffix(x) == s and prefix(x) == "rna" for x in roles.scores
                ):
                    roles.scores.append(fid)
    # prefer rna-prefixed scores first
    roles.scores.sort(key=lambda f: 0 if prefix(f) == "rna" else 1)

    if roles.rna:
        roles.gene_wrappers = resolve_gene_wrappers(roles.rna)
    if roles.protein:
        roles.protein_wrappers = resolve_protein_wrappers(roles.protein)
    return roles


def query(ds_name: str, max_items: int) -> dict[str, Any]:
    return {"linkedDsName": ds_name, "maxItems": max_items, "type": "RowsAsColsQuery"}


def design_fields(roles: Roles) -> list[str]:
    out: list[str] = []
    for fid in (
        roles.workstream,
        roles.tissue,
        roles.disease,
        roles.treatment,
        roles.inflammation,
        roles.sex,
        roles.cell_type,
    ):
        if not fid:
            continue
        col = roles.columns.get(fid) or {}
        nv = n_values(col)
        if col.get("datatype") not in CATEGORICAL_DTYPES:
            continue
        if nv is not None and nv > 40:
            continue
        if fid not in out:
            out.append(fid)
    return out[:6]


def default_gene(roles: Roles) -> tuple[str, str] | None:
    for gene in FEATURE_UMAP_GENES + MARKER_PANEL:
        if gene in roles.gene_wrappers:
            return gene, roles.gene_wrappers[gene]
    if roles.gene_wrappers:
        gene, wr = next(iter(roles.gene_wrappers.items()))
        return gene, wr
    return None


def default_protein(roles: Roles) -> tuple[str, str] | None:
    if roles.protein_wrappers:
        name, wr = next(iter(roles.protein_wrappers.items()))
        return name, wr
    return None


def view_study_overview(roles: Roles) -> dict[str, Any] | None:
    cells: list[dict[str, Any]] = []
    bits = [
        f"{roles.n_cells or '?'} cells in `{roles.obs}`.",
        f"RNA link: {roles.rna.ds_name}/{roles.rna.subgroup} ({len(roles.rna.names)} features)."
        if roles.rna
        else "No RNA expression link.",
        f"Protein link: {roles.protein.ds_name}/{roles.protein.subgroup} ({len(roles.protein.names)} features)."
        if roles.protein
        else "No protein expression link.",
        "Spatial coordinates are available; use the existing default view for the image/Viv overlay."
        if roles.embedding and roles.embedding[0] in {"x", "spatial_1", "global_1"}
        else "This project is dissociated single-cell data; spatial localisation is not available unless coordinates exist.",
    ]
    cells.append(textbox("Project summary", " ".join(bits), [0, 0], [12, 2]))
    design = design_fields(roles)
    for i, fid in enumerate(design[:4]):
        cells.append(row(suffix(fid), fid, [i * 3, 2], [3, 3]))
    y = 5
    if roles.disease and roles.tissue:
        cells.append(
            stacked(
                f"{suffix(roles.disease)} by {suffix(roles.tissue)}",
                [roles.disease, roles.tissue],
                [0, y],
                [6, 4],
            )
        )
        x = 6
    else:
        x = 0
    for i, fid in enumerate(roles.qc[:2]):
        mm = minmax(roles.columns[fid])
        if not mm:
            continue
        cells.append(histogram(suffix(fid), fid, mm[0], mm[1], [x + i * 3, y], [3, 4]))
    filt = [f for f in (roles.workstream, roles.tissue, roles.disease) if f]
    if filt:
        cells.append(selection("Filter cohort", filt, [0, 9], [6, 3]))
    extra = [f for f in design[4:] if f]
    if extra:
        cells.append(row(suffix(extra[0]), extra[0], [6, 9], [6, 3]))
    if len(cells) <= 1:
        roles.skipped.append("1 Study overview (no design/QC fields)")
        return None
    return {
        "dataSources": {roles.obs: {"layout": "gridstack", "panelWidth": 100}},
        "initialCharts": {roles.obs: cells},
    }


def view_atlas(roles: Roles) -> dict[str, Any] | None:
    if not roles.embedding:
        roles.skipped.append("2 Cell atlas (no UMAP)")
        return None
    xy = roles.embedding
    cells: list[dict[str, Any]] = []
    colors = [
        (roles.cell_type, "Cell type"),
        (roles.broad_type, "Broad type"),
        (roles.tissue, "Tissue"),
        (roles.disease, "Disease"),
    ]
    placed = 0
    for fid, title in colors:
        if not fid:
            continue
        col = placed % 2
        row_i = placed // 2
        cells.append(umap(title, fid, [col * 6, row_i * 5], [6, 5], xy))
        placed += 1
    y = ((placed + 1) // 2) * 5
    if roles.tissue and roles.cell_type:
        cells.append(
            stacked(
                f"{suffix(roles.cell_type)} by {suffix(roles.tissue)}",
                [roles.tissue, roles.cell_type],
                [0, y],
                [6, 4],
            )
        )
    disease = roles.disease
    cat = roles.broad_type or roles.cell_type
    if disease and cat:
        cells.append(
            stacked(
                f"{suffix(cat)} by {suffix(disease)}",
                [disease, cat],
                [6, y],
                [6, 4],
            )
        )
    if roles.cell_type:
        cells.append(row(suffix(roles.cell_type), roles.cell_type, [0, y + 4], [6, 3]))
    if not cells:
        roles.skipped.append("2 Cell atlas (no colour fields)")
        return None
    return {
        "dataSources": {roles.obs: {"layout": "gridstack", "panelWidth": 100}},
        "initialCharts": {roles.obs: cells},
    }


def view_markers(roles: Roles, feature_ds: dict[str, Any] | None) -> dict[str, Any] | None:
    if not roles.rna:
        roles.skipped.append("3 Marker genes (no RNA link)")
        return None
    pinned = list(roles.gene_wrappers.values())
    markers = roles.markers
    top_hits = unique_top_markers(markers)
    top_wraps: list[str] = []
    if roles.rna and top_hits:
        top_wraps = list(resolve_named_wrappers(roles.rna, [g for g, _c in top_hits]).values())
    expr_wraps = top_wraps if len(top_wraps) >= 3 else pinned
    cells: list[dict[str, Any]] = []
    y = 0
    note = (
        "Top 2 markers per cluster (unique, cap 24) on feature plots, plus the "
        f"top-20 table from `{roles.rna.matrix or roles.rna.subgroup}`. "
        "Factor dots show genes whose mean expression shifts within a cell type."
        if top_hits
        else "Canonical lineage panel. Computed top-20 / factor-varying tables were not added."
    )
    cells.append(textbox("Marker genes", note, [0, y], [12, 2]))
    y += 2
    if roles.cell_type and len(expr_wraps) >= 3:
        label = "Top cluster markers" if top_wraps else "Canonical markers"
        cells.append(
            dot(
                f"{label} by {suffix(roles.cell_type)}",
                [roles.cell_type, *expr_wraps],
                [0, y],
                [6, 5],
            )
        )
        cells.append(
            heatmap(
                f"Marker heatmap ({suffix(roles.cell_type)})",
                [roles.cell_type, *expr_wraps],
                [6, y],
                [6, 5],
            )
        )
        y += 5
    elif roles.cell_type:
        cells.append(
            dot(
                f"Gene collection by {suffix(roles.cell_type)}",
                [roles.cell_type, query(roles.rna.ds_name, 10)],
                [0, y],
                [12, 5],
            )
        )
        y += 5
    placed = 0
    if markers and roles.rna:
        for fid, genes in markers.factor_genes.items():
            wraps = list(resolve_named_wrappers(roles.rna, genes).values())
            if len(wraps) < 3:
                continue
            cells.append(
                dot(
                    f"Genes varying by {suffix(fid)}",
                    [fid, *wraps],
                    [(placed % 2) * 6, y + (placed // 2) * 5],
                    [6, 5],
                )
            )
            placed += 1
        if placed:
            y += ((placed + 1) // 2) * 5
    axes = [
        (roles.tissue, "tissue"),
        (roles.disease, "diagnosis"),
        (roles.treatment, "treatment"),
        (roles.inflammation, "inflammation"),
        (roles.response, "response"),
        (roles.workstream, "workstream"),
    ]
    if len(expr_wraps) >= 3:
        extra_dots = 0
        used = set(markers.factor_genes) if markers else set()
        axis_label = "Top cluster markers" if top_wraps else "Canonical markers"
        for fid, label in axes:
            if not fid or fid in used:
                continue
            cells.append(
                dot(
                    f"{axis_label} by {label}",
                    [fid, *expr_wraps],
                    [(extra_dots % 2) * 6, y + (extra_dots // 2) * 5],
                    [6, 5],
                )
            )
            extra_dots += 1
        y += ((extra_dots + 1) // 2) * 5 if extra_dots else 0
    if roles.embedding:
        shown = 0
        feature_genes: list[tuple[str, str]] = list(top_hits)
        if not feature_genes:
            feature_genes = [(g, "") for g in FEATURE_UMAP_GENES if g in roles.gene_wrappers][:3]
        wraps = resolve_named_wrappers(roles.rna, [g for g, _c in feature_genes]) if roles.rna else {}
        for gene, cluster in feature_genes:
            wr = wraps.get(gene) or roles.gene_wrappers.get(gene)
            if not wr:
                continue
            title = f"{gene} ({cluster})" if cluster else f"{gene} expression"
            cells.append(
                umap(
                    title,
                    wr,
                    [(shown % 4) * 3, y + (shown // 4) * 4],
                    [3, 4],
                    roles.embedding,
                )
            )
            shown += 1
        y += ((shown + 3) // 4) * 4 if shown else 0
    if roles.embedding:
        for i, fid in enumerate(roles.scores[:3]):
            cells.append(umap(suffix(fid), fid, [i * 4, y], [4, 4], roles.embedding))
    extra: dict[str, list[dict[str, Any]]] = {roles.obs: cells}
    widths = {roles.obs: {"layout": "gridstack", "panelWidth": 100}}
    if markers and markers.rows and markers.table_fields:
        tcols = [c for c in markers.table_fields if c]
        marker_charts = [
            selection("Filter markers", ["contrast", "cluster"], [0, 0], [12, 4]),
            table("Top 20 markers / factor-varying genes", tcols, [0, 4], [12, 8]),
        ]
        extra[MARKER_DS_NAME] = marker_charts
        widths[roles.obs] = {"layout": "gridstack", "panelWidth": 70}
        widths[MARKER_DS_NAME] = {"layout": "gridstack", "panelWidth": 30}
    elif feature_ds:
        avail = {c.get("field") for c in feature_ds.get("columns") or []}
        tcols = [c for c in ("name", "mean", "std", "highly_variable", "mean_counts") if c in avail]
        genes_charts: list[dict[str, Any]] = []
        if "name" in avail:
            genes_charts.append(selection("Search features", ["name"], [0, 0], [12, 4]))
        if tcols:
            genes_charts.append(table("Feature table", tcols, [0, 4], [12, 8]))
        if genes_charts:
            extra[roles.rna.ds_name] = genes_charts
            widths[roles.obs] = {"layout": "gridstack", "panelWidth": 70}
            widths[roles.rna.ds_name] = {"layout": "gridstack", "panelWidth": 30}
    if len(cells) <= 1:
        roles.skipped.append("3 Marker genes (no cell-type or markers)")
        return None
    return {"dataSources": widths, "initialCharts": extra}


def view_condition(roles: Roles) -> dict[str, Any] | None:
    if not (roles.disease or roles.tissue):
        roles.skipped.append("4 Tissue and disease (no disease/tissue)")
        return None
    cells: list[dict[str, Any]] = []
    if roles.embedding:
        if roles.disease:
            cells.append(umap("Disease", roles.disease, [0, 0], [4, 5], roles.embedding))
        if roles.inflammation:
            cells.append(umap("Inflammation", roles.inflammation, [4, 0], [4, 5], roles.embedding))
        elif roles.treatment:
            cells.append(umap("Treatment", roles.treatment, [4, 0], [4, 5], roles.embedding))
        if roles.tissue:
            cells.append(umap("Tissue", roles.tissue, [8, 0], [4, 5], roles.embedding))
    if roles.cell_type and roles.disease:
        cells.append(
            stacked(
                f"{suffix(roles.cell_type)} by {suffix(roles.disease)}",
                [roles.disease, roles.cell_type],
                [0, 5],
                [6, 4],
            )
        )
    if roles.cell_type and roles.tissue:
        cells.append(
            stacked(
                f"{suffix(roles.cell_type)} by {suffix(roles.tissue)}",
                [roles.tissue, roles.cell_type],
                [6, 5],
                [6, 4],
            )
        )
    if roles.sample_id and roles.tissue and roles.cell_type:
        cells.append(
            abundance(
                f"Per-sample {suffix(roles.cell_type)} by {suffix(roles.tissue)}",
                [roles.tissue, roles.sample_id, roles.cell_type],
                [0, 9],
                [6, 4],
            )
        )
    gene = default_gene(roles)
    if gene and roles.disease:
        cells.append(violin(f"{gene[0]} by {suffix(roles.disease)}", [roles.disease, gene[1]], [0, 13], [4, 4]))
    if gene and roles.tissue:
        cells.append(violin(f"{gene[0]} by {suffix(roles.tissue)}", [roles.tissue, gene[1]], [4, 13], [4, 4]))
    if roles.rna and roles.disease:
        cells.append(
            dot(
                f"Gene collection by {suffix(roles.disease)}",
                [roles.disease, query(roles.rna.ds_name, 10)],
                [8, 13],
                [4, 4],
            )
        )
    filt = [
        f
        for f in (
            roles.disease,
            roles.tissue,
            roles.workstream,
            roles.inflammation,
            roles.treatment,
            roles.response,
            roles.sex,
            roles.cell_type,
        )
        if f
    ]
    if filt:
        cells.append(selection("Filter disease / tissue / group", filt, [0, 17], [12, 3]))
    if not cells:
        roles.skipped.append("4 Tissue and disease (nothing to plot)")
        return None
    return {
        "dataSources": {roles.obs: {"layout": "gridstack", "panelWidth": 100}},
        "initialCharts": {roles.obs: cells},
    }


def view_rna_protein(roles: Roles, rna_ds: dict[str, Any] | None, prot_ds: dict[str, Any] | None) -> dict[str, Any] | None:
    if not roles.rna or not roles.protein:
        roles.skipped.append("5 RNA-protein concordance (no protein link)")
        return None
    gene = default_gene(roles)
    prot = default_protein(roles)
    gene_label = gene[0] if gene else "RNA"
    prot_label = prot[0] if prot else "protein"
    cells = [
        textbox(
            "RNA vs protein",
            f"Compare `{roles.rna.subgroup}` on `{roles.rna.ds_name}` with `{roles.protein.subgroup}` on `{roles.protein.ds_name}`. "
            f"UMAPs are coloured by {gene_label} RNA and {prot_label} protein when wrappers resolved. "
            "Repertoire/empty links are omitted.",
            [0, 0],
            [12, 2],
        )
    ]
    cat = roles.cell_type
    if cat:
        cells.append(
            dot(f"RNA by {suffix(cat)}", [cat, query(roles.rna.ds_name, 10)], [0, 2], [6, 4])
        )
        cells.append(
            dot(f"Protein by {suffix(cat)}", [cat, query(roles.protein.ds_name, 10)], [6, 2], [6, 4])
        )
    if roles.embedding and gene:
        cells.append(umap(f"RNA UMAP ({gene_label} RNA)", gene[1], [0, 6], [4, 4], roles.embedding))
    if roles.embedding and prot:
        cells.append(umap(f"RNA UMAP ({prot_label} protein)", prot[1], [4, 6], [4, 4], roles.embedding))
    if roles.protein_embedding and prot:
        cells.append(
            umap(
                f"Protein UMAP ({prot_label})",
                prot[1],
                [8, 6],
                [4, 4],
                roles.protein_embedding,
            )
        )
    extra: dict[str, list[dict[str, Any]]] = {roles.obs: cells}
    widths = {roles.obs: {"layout": "gridstack", "panelWidth": 70}}
    if rna_ds:
        avail = {c.get("field") for c in rna_ds.get("columns") or []}
        tcols = [c for c in ("name", "mean", "highly_variable", "mean_counts") if c in avail]
        charts = [selection("Search genes", ["name"], [0, 0], [12, 4])] if "name" in avail else []
        if tcols:
            charts.append(table("Genes", tcols, [0, 4], [12, 8]))
        if charts:
            extra[roles.rna.ds_name] = charts
            widths[roles.rna.ds_name] = {"layout": "gridstack", "panelWidth": 15}
    if prot_ds:
        avail = {c.get("field") for c in prot_ds.get("columns") or []}
        tcols = [c for c in ("name", "Protein", "gene_symbols", "highly_variable") if c in avail]
        charts = [selection("Search proteins", ["name"], [0, 0], [12, 4])] if "name" in avail else []
        if tcols:
            charts.append(table("Proteins", tcols, [0, 4], [12, 8]))
        if charts:
            extra[roles.protein.ds_name] = charts
            widths[roles.protein.ds_name] = {"layout": "gridstack", "panelWidth": 15}
    return {"dataSources": widths, "initialCharts": extra}


def print_inventory(roles: Roles) -> None:
    print("Inferred roles:")
    print(f"  obs={roles.obs} n_cells={roles.n_cells}")
    if roles.rna:
        print(
            f"  rna={roles.rna.ds_name}/{roles.rna.subgroup} n={len(roles.rna.names)} "
            f"matrix={roles.rna.matrix}"
        )
    if roles.protein:
        print(f"  protein={roles.protein.ds_name}/{roles.protein.subgroup} n={len(roles.protein.names)}")
    for label, val in [
        ("cell_type", roles.cell_type),
        ("broad_type", roles.broad_type),
        ("tissue", roles.tissue),
        ("disease", roles.disease),
        ("treatment", roles.treatment),
        ("inflammation", roles.inflammation),
        ("workstream", roles.workstream),
        ("sex", roles.sex),
        ("sample_id", roles.sample_id),
        ("embedding", roles.embedding),
        ("protein_embedding", roles.protein_embedding),
        ("qc", roles.qc),
        ("scores", roles.scores[:6]),
        ("factors", roles.factors),
    ]:
        print(f"  {label}={val}")
    print(f"  pinned_genes={list(roles.gene_wrappers)[:12]}")
    print(f"  pinned_proteins={list(roles.protein_wrappers)[:8]}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Add inferred default views to an MDV project.")
    parser.add_argument("--project", required=True, help="Path to an existing MDV project directory.")
    parser.add_argument(
        "--skip-markers",
        action="store_true",
        help="Do not compute or attach top-20 / factor-varying marker tables.",
    )
    parser.add_argument(
        "--recompute-markers",
        action="store_true",
        help="Recompute cluster and factor markers even if a cache exists.",
    )
    return parser.parse_args()


def create_default_views(
    project_path: Path | str,
    *,
    skip_markers: bool = False,
    recompute_markers: bool = False,
) -> dict[str, Any]:
    """Infer roles and write the five template views into an MDV project.

    Returns a summary with added view names and skip reasons.
    """
    project_path = Path(project_path).expanduser().resolve()
    views_path = project_path / "views.json"
    state_path = project_path / "state.json"
    ds_path = project_path / "datasources.json"
    if not views_path.is_file() or not state_path.is_file() or not ds_path.is_file():
        raise FileNotFoundError(f"Not a complete MDV project: {project_path}")

    bak = project_path / "views.json.bak"
    if not bak.exists():
        shutil.copy2(views_path, bak)
        print(f"Backed up views.json to {bak}")
    else:
        print(f"Left existing backup in place: {bak}")

    datasources = load_json(ds_path)
    roles = infer_roles(datasources)
    print_inventory(roles)
    roles.markers = load_or_compute_markers(
        project_path, roles, skip=skip_markers, recompute=recompute_markers
    )

    datasources = load_json(ds_path)
    rna_ds = ds_by_name(datasources, roles.rna.ds_name) if roles.rna else None
    prot_ds = ds_by_name(datasources, roles.protein.ds_name) if roles.protein else None

    builders = {
        "1 Study overview": lambda: view_study_overview(roles),
        "2 Cell atlas": lambda: view_atlas(roles),
        "3 Marker genes and signatures": lambda: view_markers(roles, rna_ds),
        "4 Tissue and disease context": lambda: view_condition(roles),
        "5 RNA-protein concordance": lambda: view_rna_protein(roles, rna_ds, prot_ds),
    }
    new_views: dict[str, Any] = {}
    for name, build in builders.items():
        view = build()
        if view:
            new_views[name] = view

    views = load_json(views_path)
    if not isinstance(views, dict):
        raise TypeError("views.json must be an object")
    views.update(new_views)
    save_json(views_path, views)

    state = load_json(state_path)
    all_views = list(state.get("all_views") or [])
    for name in VIEW_NAMES:
        if name in new_views and name not in all_views:
            all_views.append(name)
    state["all_views"] = all_views
    save_json(state_path, state)

    print(f"Wrote {len(new_views)} views into {project_path}")
    for name in VIEW_NAMES:
        mark = "added" if name in new_views else "skipped"
        print(f"  - {name} ({mark})")
    if roles.skipped:
        print("Skip reasons:")
        for s in roles.skipped:
            print(f"  - {s}")
    return {"added": list(new_views), "skipped": list(roles.skipped)}


def main() -> None:
    args = parse_args()
    create_default_views(
        args.project,
        skip_markers=args.skip_markers,
        recompute_markers=args.recompute_markers,
    )


if __name__ == "__main__":
    main()
