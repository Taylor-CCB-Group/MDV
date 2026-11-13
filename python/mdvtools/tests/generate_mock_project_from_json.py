#!/usr/bin/env python3
"""Build a mock MDV project from a project's JSON files alone.

Copies views.json and state.json unchanged and writes a datafile.h5 with made-up values for
every column in datasources.json, so views can be checked locally without the real h5.

- Cells (any table that is not the feature side of a rows-as-columns link) are capped at
  --max-rows; only that table's ``size`` changes in the mock datasources.json.
- Feature tables (genes, proteins) keep their size and row order, so gene positions saved in
  views (``subgroup|GENE (subgroup)|row``) still point at the right gene.
- Rows share a hidden cluster, so categories and numeric columns (UMAPs included) form blobs,
  and genes named in the views get cluster-specific expression.

Values are random. Use it to check layout, links and legends, not biology.

Usage:
    python -m mdvtools.tests.generate_mock_project_from_json \\
        --source ~/Documents/MDV_resources/synovial_latest --out ~/mdv/mock_synovial_latest
Then open http://localhost:5055/rescan_projects (or restart the server) to register it.
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
from pathlib import Path
from typing import Any

import h5py
import numpy as np

MARKER = ".mock_project"
CLUSTERS = 24
NOISE = 0.15
FIXED_FEATURE = re.compile(r"^([^|]+)\|(.+) \(\1\)\|(\d+)$")


def load_json(path: Path) -> Any:
    with path.open() as handle:
        return json.load(handle)


def rows_as_columns_links(ds: dict[str, Any]) -> dict[str, dict[str, Any]]:
    """Linked datasource name -> rows_as_columns spec, for links this datasource owns."""
    links = ds.get("links") or {}
    return {
        target: link["rows_as_columns"]
        for target, link in links.items()
        if isinstance(link, dict) and isinstance(link.get("rows_as_columns"), dict)
    }


def walk_strings(value: Any):
    if isinstance(value, str):
        yield value
    elif isinstance(value, dict):
        for item in value.values():
            yield from walk_strings(item)
    elif isinstance(value, list):
        for item in value:
            yield from walk_strings(item)


def features_named_in_views(views: dict[str, Any], feature_ds: set[str], subgroup_ds: dict[str, str]):
    """Return ({feature_ds: {row: name}} from fixed fields, {feature_ds: {names}} worth expressing)."""
    fixed: dict[str, dict[int, str]] = {name: {} for name in feature_ds}
    wanted: dict[str, set[str]] = {name: set() for name in feature_ds}
    for view in views.values():
        if not isinstance(view, dict):
            continue
        charts_by_ds = view.get("initialCharts") or {}
        for text in walk_strings(charts_by_ds):
            match = FIXED_FEATURE.match(text)
            if match and match.group(1) in subgroup_ds:
                target = subgroup_ds[match.group(1)]
                fixed[target][int(match.group(3))] = match.group(2)
                wanted[target].add(match.group(2))
        for ds_name, charts in charts_by_ds.items():
            if ds_name not in feature_ds:
                continue
            for chart in charts:
                for spec in (chart.get("filters") or {}).values():
                    if isinstance(spec, dict):
                        wanted[ds_name].update(v for v in spec.get("category") or [] if isinstance(v, str))
    return fixed, wanted


class MockWriter:
    def __init__(self, seed: int, background: int):
        self.rng = np.random.default_rng(seed)
        self.background = background

    def clusters(self, n: int) -> np.ndarray:
        weights = self.rng.dirichlet(np.full(CLUSTERS, 0.8))
        return self.rng.choice(CLUSTERS, size=n, p=weights)

    def noisy(self, mapped: np.ndarray, high: int) -> np.ndarray:
        flip = self.rng.random(mapped.size) < NOISE
        mapped = mapped.copy()
        mapped[flip] = self.rng.integers(0, high, flip.sum())
        return mapped

    def column(self, col: dict[str, Any], n: int, cl: np.ndarray, one_per_row: bool):
        """Return (data, dtype) for one column, or None for a datatype this script can't fill."""
        dt = col.get("datatype")
        if dt in ("text", "text16"):
            values = col.get("values") or ["ND"]
            col["values"] = values
            high = len(values)
            if one_per_row and high >= n:
                raw = np.arange(n)
            else:
                raw = self.noisy(self.rng.integers(0, high, CLUSTERS)[cl], high)
            return raw, (np.uint8 if dt == "text" and high <= 256 else np.uint16)
        if dt == "multitext":
            values = col.get("values") or ["ND"]
            col["values"] = values
            width = max(int(col.get("stringLength") or 1), 1)
            col["stringLength"] = width
            raw = np.full((n, width), 65535, dtype=np.uint16)
            raw[:, 0] = self.noisy(self.rng.integers(0, len(values), CLUSTERS)[cl], len(values))
            return raw.ravel(), np.uint16
        if dt == "unique":
            strings = [f"{col['field']}_{i}" for i in range(n)]
            width = max(int(col.get("stringLength") or 0), max(len(s) for s in strings))
            col["stringLength"] = width
            return strings, h5py.string_dtype("utf-8", width)
        if dt in ("double", "integer", "int32"):
            low, high = (col.get("minMax") or [0.0, 1.0])[:2]
            if not high > low:
                high = low + 1.0
            centers = self.rng.uniform(low, high, CLUSTERS)
            data = np.clip(centers[cl] + self.rng.normal(0, (high - low) * 0.06, n), low, high)
            if dt == "int32":
                return np.rint(data), np.int32
            if dt == "integer":
                data = np.rint(data)
            return data, np.float32
        return None

    def expression(self, n: int, cl: np.ndarray, n_features: int, rich: set[int]):
        """Per-feature (cell indices, log-scale values); rich features vary by cluster."""
        k = min(self.background, n)
        background = np.sort(self.rng.integers(0, n, (n_features, k)), axis=1).astype(np.uint32)
        out: list[tuple[np.ndarray, np.ndarray]] = []
        for feature in range(n_features):
            if feature in rich:
                p_on = self.rng.beta(0.6, 1.4, CLUSTERS)
                level = self.rng.uniform(0.5, 3.0, CLUSTERS)
                idx = np.flatnonzero(self.rng.random(n) < p_on[cl]).astype(np.uint32)
                vals = np.clip(level[cl[idx]] + self.rng.normal(0, 0.5, idx.size), 0.05, None)
            else:
                idx = np.unique(background[feature])
                vals = self.rng.uniform(0.1, 1.5, idx.size)
            out.append((idx, vals.astype(np.float32)))
        return out


def write_sparse(group: h5py.Group, per_feature: list[tuple[np.ndarray, np.ndarray]], counts: bool):
    sizes = np.fromiter((idx.size for idx, _ in per_feature), dtype=np.int64, count=len(per_feature))
    pointers = np.concatenate([[0], np.cumsum(sizes)])
    indices = np.concatenate([idx for idx, _ in per_feature]) if per_feature else np.array([], np.uint32)
    values = np.concatenate([vals for _, vals in per_feature]) if per_feature else np.array([], np.float32)
    if counts:
        values = np.maximum(1, np.rint(np.expm1(values))).astype(np.float32)
    group.create_dataset("x", data=values, dtype=np.float32)
    group.create_dataset("i", data=indices, dtype=np.uint32)
    group.create_dataset("p", data=pointers)


def write_dense(group: h5py.Group, n: int, per_feature: list[tuple[np.ndarray, np.ndarray]], counts: bool):
    dense = np.zeros(n * len(per_feature), dtype=np.float32)
    for feature, (idx, vals) in enumerate(per_feature):
        dense[feature * n + idx.astype(np.int64)] = np.maximum(1, np.rint(np.expm1(vals))) if counts else vals
    group.create_dataset("x", data=dense)
    group.create_dataset("length", data=[n])


def build(source: Path, out: Path, max_rows: int, seed: int, background: int) -> dict[str, Any]:
    datasources = load_json(source / "datasources.json")
    if isinstance(datasources, dict):
        datasources = [datasources]
    views = load_json(source / "views.json")
    state = load_json(source / "state.json")

    by_name = {ds["name"]: ds for ds in datasources}
    feature_ds: set[str] = set()
    subgroup_ds: dict[str, str] = {}
    for ds in datasources:
        for target, spec in rows_as_columns_links(ds).items():
            feature_ds.add(target)
            for stub in spec.get("subgroups") or {}:
                subgroup_ds[stub] = target
    fixed, wanted = features_named_in_views(views, feature_ds, subgroup_ds)

    writer = MockWriter(seed, background)
    clusters: dict[str, np.ndarray] = {}
    name_rows: dict[str, dict[str, int]] = {}
    skipped: list[str] = []
    resized: dict[str, tuple[int, int]] = {}
    groups: dict[str, h5py.Group] = {}

    out.mkdir(parents=True, exist_ok=True)
    (out / MARKER).write_text(f"Mock project generated from {source.name}. Safe to delete.\n")
    with h5py.File(out / "datafile.h5", "w") as h5:
        for ds in datasources:
            name = ds["name"]
            n = int(ds["size"])
            if name not in feature_ds and n > max_rows:
                resized[name] = (n, max_rows)
                n = max_rows
                ds["size"] = n
            cl = writer.clusters(n)
            clusters[name] = cl
            group = h5.create_group(name)
            groups[name] = group
            name_columns = {
                spec.get("name_column", "name")
                for other in datasources
                for target, spec in rows_as_columns_links(other).items()
                if target == name
            }
            for col in ds["columns"]:
                made = writer.column(col, n, cl, one_per_row=col["field"] in name_columns)
                if made is None:
                    skipped.append(f"{name}:{col['field']} ({col.get('datatype')})")
                    continue
                data, dtype = made
                if col["field"] in name_columns and col.get("values"):
                    data = np.asarray(data)
                    lookup = {value: i for i, value in enumerate(col["values"])}
                    for row, feature in fixed.get(name, {}).items():
                        if row < n and feature in lookup:
                            other = np.flatnonzero(data == lookup[feature])
                            if other.size:
                                data[other[0]], data[row] = data[row], lookup[feature]
                    name_rows[name] = {col["values"][v]: row for row, v in enumerate(data)}
                group.create_dataset(col["field"], data=data, dtype=dtype)

        for ds in datasources:
            n = int(ds["size"])
            for target, spec in rows_as_columns_links(ds).items():
                if target not in by_name:
                    continue
                n_features = int(by_name[target]["size"])
                rows = name_rows.get(target, {})
                rich = {rows[g] for g in wanted.get(target, set()) if g in rows}
                per_feature = writer.expression(n, clusters[ds["name"]], n_features, rich)
                for stub, info in (spec.get("subgroups") or {}).items():
                    label = f"{stub} {info.get('name', '')}".lower()
                    counts = ("raw" in label or "count" in label) and "log" not in label
                    group = groups[ds["name"]].create_group(str(info.get("name", stub)))
                    if info.get("type") == "sparse" or n * n_features > 50_000_000:
                        info["type"] = "sparse"
                        write_sparse(group, per_feature, counts)
                    else:
                        write_dense(group, n, per_feature, counts)

    with (out / "datasources.json").open("w") as handle:
        json.dump(datasources, handle)
    shutil.copyfile(source / "views.json", out / "views.json")
    with (out / "state.json").open("w") as handle:
        json.dump(state, handle, indent=2)
    for path in out.iterdir():
        path.chmod(0o644)
    return {"resized": resized, "skipped": skipped, "features_named_in_views": {k: len(v) for k, v in wanted.items()}}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--source", type=Path, required=True, help="Folder with datasources.json, views.json, state.json.")
    parser.add_argument("--out", type=Path, required=True, help="New project folder, e.g. ~/mdv/mock_myeloid.")
    parser.add_argument("--max-rows", type=int, default=20_000, help="Cap for cell tables (default 20000).")
    parser.add_argument("--background", type=int, default=200, help="Cells expressing each unnamed gene (default 200).")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--force", action="store_true", help="Replace --out if it is an earlier mock project.")
    args = parser.parse_args()

    source = args.source.expanduser().resolve()
    out = args.out.expanduser().resolve()
    if out == source:
        parser.error("--out must differ from --source")
    if out.exists():
        if not (out / MARKER).exists():
            parser.error(f"{out} exists and is not a mock project; refusing to touch it")
        if not args.force:
            parser.error(f"{out} exists; pass --force to rebuild it")
        shutil.rmtree(out)
    summary = build(source, out, args.max_rows, args.seed, args.background)
    print(f"Wrote mock project to {out}")
    for name, (before, after) in summary["resized"].items():
        print(f"  {name}: {before} rows -> {after}")
    for name, count in summary["features_named_in_views"].items():
        print(f"  {name}: {count} features named in views get cluster-specific values")
    if summary["skipped"]:
        print("  Skipped columns:", ", ".join(summary["skipped"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
