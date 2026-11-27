"""Export an MDV project to one SpatialData store per image region.

Each store has the same element names: ``images/morphology``,
``shapes/cell_boundaries``, and ``tables/table``. Marker column groups become
the table's variables and layers. Image channels are pyramided until the
longest side of the coarsest level is at most 1024 px.

``export_converted_project`` writes those stores into a new MDV project whose
marker columns are a ``markers`` datasource, with views rewritten to
rows-as-columns fields.
"""

from __future__ import annotations

import json
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd

from mdvtools.mdvproject import MDVProject

IMAGE_NAME = "morphology"
SHAPES_NAME = "cell_boundaries"
TABLE_NAME = "table"
REGION_KEY = "region"
INSTANCE_KEY = "cell_id"
OBJECT_ID_FIELD = "Object Id"
SEGMENTATION_PARTS = "segmentation_parts"
PRIMARY_LAYER_NAME = "log2"
MARKERS_DATASOURCE = "markers"
MARKERS_NAME_COLUMN = "name"
_NUMERIC_SUFFIX = re.compile(r"^[0-9]+$")
_WHOLE_NUMBER_TOLERANCE = 1e-3


@dataclass(frozen=True)
class ColumnLayout:
    """How datasource columns land in one AnnData table."""

    var_names: list[str]
    layers: dict[str, list[str]]
    primary_layer: str
    obsm: dict[str, list[str]]
    obs_columns: list[str]


def parse_channel_renames(pairs: Sequence[str]) -> dict[str, str]:
    """Parse ``OLD=NEW`` channel rename pairs. Later pairs do not override earlier ones."""
    renames: dict[str, str] = {}
    for pair in pairs:
        if pair.count("=") != 1:
            raise ValueError(f"Channel rename must be OLD=NEW, got {pair!r}")
        old, new = (part.strip() for part in pair.split("=", 1))
        if not old or not new:
            raise ValueError(f"Channel rename must be OLD=NEW, got {pair!r}")
        if old in renames:
            raise ValueError(f"Duplicate channel rename for {old!r}")
        renames[old] = new
    return renames


def export_mdv_to_spatialdata(
    project_dir: str | Path,
    output_dir: str | Path,
    *,
    datasource: str = "cells",
    channel_renames: Mapping[str, str] | None = None,
) -> list[Path]:
    """Write one SpatialData zarr per region in ``datasource``.

    ``channel_renames`` maps image channel names (for example ``DAPI_FINAL`` to
    ``DAPI``) before the image is written. Table variable names are unchanged.
    A relative region ``json`` path is joined to the project directory. A path
    that does not exist raises ``FileNotFoundError``.
    """
    project_path = Path(project_dir)
    destination_root = Path(output_dir)
    renames = dict(channel_renames or {})
    project = MDVProject(str(project_path))
    metadata = project.get_datasource_metadata(datasource)
    regions = metadata.get("regions")
    if not isinstance(regions, dict) or not regions.get("all_regions"):
        raise ValueError(f"Datasource {datasource!r} has no regions")
    position_fields = [str(field) for field in regions["position_fields"]]
    region_field = str(regions["region_field"])
    project_scale = float(regions.get("scale", 1.0))
    frame = project.get_datasource_as_dataframe(datasource)
    layout = classify_columns(
        frame,
        list(metadata.get("columnGroups") or []),
        position_fields,
    )

    destination_root.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    for region_name, region in regions["all_regions"].items():
        _check_region_name(str(region_name))
        destination = destination_root / f"{region_name}.zarr"
        if destination.exists():
            raise FileExistsError(f"SpatialData store already exists: {destination}")
        _write_region(
            project_path,
            destination,
            frame,
            layout,
            regions,
            str(region_name),
            region,
            region_field=region_field,
            project_scale=project_scale,
            channel_renames=renames,
        )
        written.append(destination)
        print(f"Wrote {destination}", flush=True)
    return written


def export_converted_project(
    project_dir: str | Path,
    output_dir: str | Path,
    *,
    datasource: str = "cells",
    markers_datasource: str = MARKERS_DATASOURCE,
    channel_renames: Mapping[str, str] | None = None,
) -> Path:
    """Write an MDV project whose marker columns are a ``markers`` datasource.

    SpatialData stores are written under ``output_dir/spatial``, one zarr per
    region. The cell table keeps its other columns. Each marker layer is a
    rows-as-columns subgroup, and view fields that named a marker column
    (``log2_CD3``) become the link field (``log2|CD3 (log2)|index``).
    ``images`` and ``json`` are symlinked from the source project.
    """
    source_path = Path(project_dir)
    output_path = Path(output_dir)
    if output_path.exists() and any(output_path.iterdir()):
        raise FileExistsError(f"Output project already exists: {output_path}")
    if markers_datasource == datasource:
        raise ValueError("The markers datasource must be a different name from the cell datasource")

    source = MDVProject(str(source_path))
    metadata = source.get_datasource_metadata(datasource)
    regions = metadata.get("regions")
    if not isinstance(regions, dict) or not regions.get("all_regions"):
        raise ValueError(f"Datasource {datasource!r} has no regions")
    position_fields = [str(field) for field in regions["position_fields"]]
    frame = source.get_datasource_as_dataframe(datasource)
    column_groups = list(metadata.get("columnGroups") or [])
    layout = classify_columns(frame, column_groups, position_fields)

    project = MDVProject(str(output_path), delete_existing=False)
    kept = _columns_except_markers(frame, layout)
    project.add_datasource(datasource, frame.loc[:, kept].copy(), add_to_view=None)
    for name, columns in _groups_among(column_groups, set(kept)):
        project.set_column_group(datasource, name, columns)
    _attach_regions(project, datasource, regions)
    _symlink_project_dir(source_path, output_path, "images")
    _symlink_project_dir(source_path, output_path, "json")

    markers = pd.DataFrame({MARKERS_NAME_COLUMN: layout.var_names})
    project.add_datasource(markers_datasource, markers, add_to_view=None)
    project.add_rows_as_columns_link(datasource, markers_datasource, MARKERS_NAME_COLUMN, markers_datasource)
    for layer, columns in layout.layers.items():
        project.add_rows_as_columns_subgroup(
            datasource,
            markers_datasource,
            layer,
            _matrix(frame, columns),
            sparse=False,
        )
    _copy_views(source, project, expression_field_map(layout))
    export_mdv_to_spatialdata(
        source_path,
        output_path / "spatial",
        datasource=datasource,
        channel_renames=channel_renames,
    )
    return output_path


def expression_field(layer: str, marker: str, index: int) -> str:
    """Rows-as-columns field string, ``layer|marker (layer)|index``."""
    return f"{layer}|{marker} ({layer})|{index}"


def expression_field_map(layout: ColumnLayout) -> dict[str, str]:
    """Map each flat marker column (``log2_CD3``) to its rows-as-columns field."""
    mapping: dict[str, str] = {}
    for layer, columns in layout.layers.items():
        for index, column in enumerate(columns):
            mapping[column] = expression_field(layer, layout.var_names[index], index)
    return mapping


def rewrite_expression_fields(value: Any, mapping: Mapping[str, str]) -> Any:
    """Replace marker column names anywhere they appear as a whole string in ``value``."""
    if isinstance(value, str):
        return mapping.get(value, value)
    if isinstance(value, list):
        return [rewrite_expression_fields(item, mapping) for item in value]
    if isinstance(value, dict):
        return {key: rewrite_expression_fields(item, mapping) for key, item in value.items()}
    return value


def classify_columns(
    frame: pd.DataFrame,
    column_groups: Sequence[Mapping[str, object]],
    position_fields: Sequence[str],
) -> ColumnLayout:
    """Split columns into marker layers, embeddings, and observation columns.

    Marker groups are column groups whose columns are ``{group}_{marker}`` and
    whose markers are not plain integers. Every marker group must name the same
    set of markers. ``log2`` is the primary matrix when that group exists.
    Integer suffixes (``X_pca_1``) become ``obsm``. ``position_fields`` become
    ``obsm['spatial']``. A group named ``spatial`` is the same coordinates when
    its values match ``position_fields``, and stays in ``obs`` when they do not.
    """
    fields = [str(field) for field in position_fields]
    if len(fields) != 2:
        raise ValueError("regions.position_fields must be two columns, x then y")
    missing = [field for field in fields if field not in frame.columns]
    if missing:
        raise ValueError(f"Position fields missing from the datasource: {', '.join(missing)}")
    if OBJECT_ID_FIELD not in frame.columns:
        raise ValueError(f"Datasource has no {OBJECT_ID_FIELD!r} column")

    consumed = set(fields)
    obsm: dict[str, list[str]] = {"spatial": list(fields)}
    marker_groups: list[tuple[str, list[str], list[str]]] = []

    for group in column_groups:
        name = str(group["name"])
        raw_columns = group.get("columns")
        if isinstance(raw_columns, str) or not isinstance(raw_columns, Sequence):
            raise ValueError(f"Column group {name!r} has no columns list")
        group_columns = [str(column) for column in raw_columns]
        if not group_columns:
            continue
        if name == "spatial" or group_columns == fields:
            if _columns_match(frame, group_columns, fields):
                consumed.update(group_columns)
            continue
        suffixes = _group_suffixes(name, group_columns)
        if suffixes is None:
            continue
        if all(_NUMERIC_SUFFIX.fullmatch(suffix) for suffix in suffixes):
            ordered = [
                column
                for _suffix, column in sorted(zip(suffixes, group_columns), key=lambda item: int(item[0]))
            ]
            _claim(consumed, ordered, name)
            obsm[name] = ordered
            continue
        if any(_NUMERIC_SUFFIX.fullmatch(suffix) for suffix in suffixes):
            continue
        if len(suffixes) != len(set(suffixes)):
            raise ValueError(f"Column group {name!r} repeats a marker")
        marker_groups.append((name, suffixes, group_columns))

    if not marker_groups:
        raise ValueError("No marker column groups found")

    var_names = marker_groups[0][1]
    var_set = set(var_names)
    layers: dict[str, list[str]] = {}
    for name, suffixes, group_columns in marker_groups:
        if set(suffixes) != var_set:
            raise ValueError(
                "Marker column groups do not share a marker set: "
                f"{marker_groups[0][0]} ({', '.join(var_names)}), "
                f"{name} ({', '.join(suffixes)})"
            )
        by_marker = dict(zip(suffixes, group_columns))
        ordered = [by_marker[marker] for marker in var_names]
        _claim(consumed, ordered, name)
        layers[name] = ordered

    consumed.add(OBJECT_ID_FIELD)
    obs_columns: list[str] = []
    for column in frame.columns:
        if column in consumed:
            continue
        if column in {INSTANCE_KEY, REGION_KEY, SEGMENTATION_PARTS}:
            raise ValueError(
                f"Column {column!r} clashes with the table annotation column of the same name"
            )
        obs_columns.append(str(column))

    primary = PRIMARY_LAYER_NAME if PRIMARY_LAYER_NAME in layers else marker_groups[0][0]
    return ColumnLayout(var_names, layers, primary, obsm, obs_columns)


def _write_region(
    project_path: Path,
    destination: Path,
    frame: pd.DataFrame,
    layout: ColumnLayout,
    regions: Mapping[str, object],
    region_name: str,
    region: Mapping[str, object],
    *,
    region_field: str,
    project_scale: float,
    channel_renames: Mapping[str, str],
) -> None:
    from spatialdata import SpatialData
    from spatialdata.models import Image2DModel, ShapesModel, TableModel
    from spatialdata_js_util import ensure_pyramids

    json_path = resolve_project_file(project_path, _region_json(region_name, region))
    image_path = resolve_project_file(project_path, _region_image(regions, region_name, region))
    region_rows = frame.loc[frame[region_field].astype(str) == region_name].copy()
    if region_rows.empty:
        raise ValueError(f"Region {region_name!r} has no cells in {region_field!r}")
    cell_ids = _cell_ids(region_rows[OBJECT_ID_FIELD], region_name)
    order = np.argsort(cell_ids, kind="mergesort")
    cell_ids = cell_ids[order]
    region_rows = region_rows.iloc[order]
    polygons = _read_segmentation(json_path)
    _check_id_join(region_name, cell_ids, polygons)
    part_counts = [polygons[int(cell_id)].parts for cell_id in cell_ids]

    image = _read_morphology(image_path)
    scale_x = image.scale_x if image.scale_x is not None else project_scale
    scale_y = image.scale_y if image.scale_y is not None else project_scale
    channel_names = _rename_channels(image.channel_names, channel_renames, image_path)
    transform = _pixel_scale(scale_x, scale_y)
    parsed_image = Image2DModel.parse(
        image.array,
        dims=("c", "y", "x"),
        c_coords=channel_names,
        transformations=transform,
    )
    shapes = ShapesModel.parse(_shapes_frame(cell_ids, polygons), transformations=_pixel_scale(scale_x, scale_y))
    table = TableModel.parse(
        _region_table(region_rows, layout, cell_ids, part_counts),
        region=SHAPES_NAME,
        region_key=REGION_KEY,
        instance_key=INSTANCE_KEY,
    )
    spatial = SpatialData(
        images={IMAGE_NAME: parsed_image},
        shapes={SHAPES_NAME: shapes},
        tables={TABLE_NAME: table},
    )
    ensure_pyramids(spatial, image_keys=[IMAGE_NAME])
    # One channel plane is one chunk. A synchronous scheduler keeps a single plane in memory.
    import dask.config

    with dask.config.set(scheduler="synchronous"):
        spatial.write(destination)


def resolve_project_file(project_dir: Path, raw: str) -> Path:
    """Resolve ``raw`` against the project when it is relative, and require the file."""
    path = Path(raw)
    resolved = path if path.is_absolute() else project_dir / path
    if not resolved.is_file():
        raise FileNotFoundError(f"Project file not found: {resolved}")
    return resolved


@dataclass(frozen=True)
class _Morphology:
    array: Any
    channel_names: list[str]
    scale_x: float | None
    scale_y: float | None


def _read_morphology(path: Path) -> _Morphology:
    """Read an OME-TIFF as a lazy ``(c, y, x)`` array.

    tifffile's zarr store does not import against the zarr installed here, so
    each channel plane is a delayed read of that page.
    """
    import dask.array as da
    import tifffile

    with tifffile.TiffFile(path) as tiff:
        series = tiff.series[0]
        axes = series.axes.upper()
        sizes = dict(zip(axes, (int(size) for size in series.shape)))
        for axis in ("T", "Z"):
            if sizes.get(axis, 1) != 1:
                raise ValueError(
                    f"{path} has {axis} size {sizes[axis]}; only a single 2D plane is supported"
                )
        if "Y" not in sizes or "X" not in sizes:
            raise ValueError(f"{path} series axes {series.axes!r} have no Y and X")
        n_channels = sizes.get("C", 1)
        height = sizes["Y"]
        width = sizes["X"]
        series_pages = series.pages
        if series_pages is None:
            raise ValueError(f"{path} has no TIFF pages")
        pages = [page for page in series_pages if page is not None]
        if len(pages) != len(series_pages):
            raise ValueError(f"{path} has a missing TIFF page")
        page_shapes = [tuple(int(dim) for dim in page.shape) for page in pages]
        dtype = np.dtype(series.dtype)
        channel_names, scale_x, scale_y = _ome_image_metadata(tiff.ome_metadata, n_channels, path)

    path_str = str(path)
    if len(pages) == n_channels and all(shape == (height, width) for shape in page_shapes):
        planes = [
            da.from_delayed(
                _delayed_page(path_str, index),
                shape=(height, width),
                dtype=dtype,
            )
            for index in range(n_channels)
        ]
    elif len(pages) == 1 and page_shapes[0] == (n_channels, height, width):
        planes = [
            da.from_delayed(
                _delayed_channel(path_str, index),
                shape=(height, width),
                dtype=dtype,
            )
            for index in range(n_channels)
        ]
    else:
        raise ValueError(
            f"Cannot read {path} as (c, y, x); series axes {axes} shape {_series_shape_text(sizes)}, "
            f"{len(pages)} pages"
        )
    return _Morphology(da.stack(planes, axis=0), channel_names, scale_x, scale_y)


def _series_shape_text(sizes: Mapping[str, int]) -> str:
    return "(" + ", ".join(f"{axis}={size}" for axis, size in sizes.items()) + ")"


def _read_series_page(path: str, page_index: int) -> np.ndarray:
    import tifffile

    with tifffile.TiffFile(path) as tiff:
        pages = tiff.series[0].pages
        if pages is None:
            raise ValueError(f"{path} has no TIFF pages")
        page = pages[page_index]
        if page is None:
            raise ValueError(f"{path} is missing page {page_index}")
        return np.asarray(page.asarray())


def _delayed_page(path: str, page_index: int):
    from dask.delayed import delayed

    @delayed
    def read() -> np.ndarray:
        array = _read_series_page(path, page_index)
        if array.ndim != 2:
            raise ValueError(f"Expected a 2D plane from {path} page {page_index}, got shape {array.shape}")
        return array

    return read()


def _delayed_channel(path: str, channel_index: int):
    from dask.delayed import delayed

    @delayed
    def read() -> np.ndarray:
        array = _read_series_page(path, 0)
        if array.ndim != 3:
            raise ValueError(f"Expected a channel stack from {path}, got shape {array.shape}")
        return array[channel_index]

    return read()


def _ome_image_metadata(
    ome_xml: str | None,
    n_channels: int,
    path: Path,
) -> tuple[list[str], float | None, float | None]:
    if not ome_xml:
        return [f"channel_{index}" for index in range(n_channels)], None, None
    import xml.etree.ElementTree as ET

    root = ET.fromstring(ome_xml)
    namespace = root.tag[root.tag.find("{") + 1 : root.tag.find("}")] if root.tag.startswith("{") else ""
    prefix = f"{{{namespace}}}" if namespace else ""
    pixels = root.find(f".//{prefix}Pixels")
    channels = [] if pixels is None else list(pixels.findall(f"{prefix}Channel"))
    names = [channel.attrib.get("Name", f"channel_{index}") for index, channel in enumerate(channels)]
    if len(names) != n_channels:
        raise ValueError(
            f"{path} OME metadata lists {len(names)} channels and the image has {n_channels}"
        )
    scale_x = _optional_float(None if pixels is None else pixels.attrib.get("PhysicalSizeX"))
    scale_y = _optional_float(None if pixels is None else pixels.attrib.get("PhysicalSizeY"))
    if scale_x is None:
        scale_x = scale_y
    elif scale_y is None:
        scale_y = scale_x
    return names, scale_x, scale_y


def _optional_float(value: str | None) -> float | None:
    if value is None or value == "":
        return None
    return float(value)


def _rename_channels(names: list[str], renames: Mapping[str, str], path: Path) -> list[str]:
    missing = [old for old in renames if old not in names]
    if missing:
        raise ValueError(
            f"{path} has no channel {missing[0]!r} to rename. Channels: {', '.join(names)}"
        )
    renamed = [renames.get(name, name) for name in names]
    if len(set(renamed)) != len(renamed):
        raise ValueError(f"Channel names are not unique after renaming {path}: {renamed}")
    return renamed


def _pixel_scale(scale_x: float, scale_y: float) -> dict[str, Any]:
    from spatialdata.transformations import Scale

    if scale_x <= 0 or scale_y <= 0:
        raise ValueError(f"Pixel scale must be positive, got x={scale_x}, y={scale_y}")
    # Declared as Any values so the dict matches Image2DModel's transformation mapping.
    transforms: dict[str, Any] = {"global": Scale([scale_x, scale_y], axes=("x", "y"))}
    return transforms


def _shapes_frame(cell_ids: np.ndarray, polygons: Mapping[int, _CellBoundary]):
    import geopandas as gpd

    instance_ids = [str(int(cell_id)) for cell_id in cell_ids]
    geometries = [polygons[int(cell_id)].geometry for cell_id in cell_ids]
    return gpd.GeoDataFrame({"geometry": geometries}, index=pd.Index(instance_ids))


def _region_table(
    frame: pd.DataFrame,
    layout: ColumnLayout,
    cell_ids: np.ndarray,
    part_counts: Sequence[int],
):
    from anndata import AnnData

    instance_ids = [str(int(cell_id)) for cell_id in cell_ids]
    obs = frame.loc[:, layout.obs_columns].copy()
    for column in obs.columns:
        if obs[column].dtype == object:
            obs[column] = obs[column].astype(str)
    obs.insert(0, REGION_KEY, pd.Categorical([SHAPES_NAME] * len(obs)))
    obs.insert(0, INSTANCE_KEY, instance_ids)
    obs.insert(2, SEGMENTATION_PARTS, np.asarray(part_counts, dtype=np.int32))
    obs.index = pd.Index(instance_ids)
    table = AnnData(
        X=_matrix(frame, layout.layers[layout.primary_layer]),
        obs=obs,
        var=pd.DataFrame(index=pd.Index(layout.var_names, dtype=str)),
        layers={name: _matrix(frame, columns) for name, columns in layout.layers.items()},
    )
    for name, columns in layout.obsm.items():
        table.obsm[name] = _matrix(frame, columns)
    return table


def _matrix(frame: pd.DataFrame, columns: Sequence[str]) -> np.ndarray:
    array = frame.loc[:, list(columns)].to_numpy(dtype=np.float32, copy=True)
    if array.ndim == 1:
        return array.reshape(-1, 1)
    return array


@dataclass(frozen=True)
class _CellBoundary:
    geometry: object
    parts: int


def _read_segmentation(path: Path) -> dict[int, _CellBoundary]:
    """Read cell polygons. Features that share a DN become one MultiPolygon."""
    from shapely.geometry import MultiPolygon, Polygon, shape

    payload = json.loads(path.read_text())
    features = payload.get("features") if isinstance(payload, dict) else None
    if not isinstance(features, list):
        raise ValueError(f"{path} has no GeoJSON features")
    parts: dict[int, list[Polygon]] = {}
    for feature in features:
        if not isinstance(feature, dict):
            raise ValueError(f"{path} has a feature that is not an object")
        properties = feature.get("properties")
        if not isinstance(properties, dict) or "DN" not in properties:
            raise ValueError(f"{path} has a feature without DN")
        dn = _whole_number(properties["DN"], f"DN in {path.name}")
        raw_geometry = feature.get("geometry")
        if not isinstance(raw_geometry, dict):
            raise ValueError(f"{path} feature DN {dn} has no geometry")
        geometry = shape(raw_geometry)
        if geometry.is_empty:
            raise ValueError(f"{path} has an empty geometry for DN {dn}")
        parts.setdefault(dn, []).extend(_polygon_parts(geometry, path, dn))
    polygons: dict[int, _CellBoundary] = {}
    for dn, pieces in parts.items():
        if len(pieces) == 1:
            polygons[dn] = _CellBoundary(pieces[0], 1)
            continue
        merged = MultiPolygon(pieces)
        if merged.is_empty or not merged.is_valid:
            raise ValueError(
                f"{path} cannot combine {len(pieces)} features with DN {dn} into a MultiPolygon"
            )
        polygons[dn] = _CellBoundary(merged, len(pieces))
    return polygons


def _polygon_parts(geometry: object, path: Path, dn: int) -> list[Any]:
    from shapely.geometry import MultiPolygon, Polygon

    if isinstance(geometry, Polygon):
        return [geometry]
    if isinstance(geometry, MultiPolygon):
        return list(geometry.geoms)
    kind = getattr(geometry, "geom_type", type(geometry).__name__)
    raise ValueError(f"{path} feature DN {dn} has geometry {kind}, expected Polygon")


def _check_id_join(region_name: str, cell_ids: np.ndarray, polygons: Mapping[int, _CellBoundary]) -> None:
    table_ids = {int(cell_id) for cell_id in cell_ids}
    segmentation_ids = set(polygons)
    if table_ids == segmentation_ids:
        return
    raise ValueError(
        f"Region {region_name!r}: GeoJSON DN and Object Id differ "
        f"({len(table_ids - segmentation_ids)} cells only in the table, "
        f"{len(segmentation_ids - table_ids)} only in the GeoJSON)"
    )


def _cell_ids(values: pd.Series, region_name: str) -> np.ndarray:
    ids = np.array(
        [_whole_number(value, f"Object Id in region {region_name!r}") for value in values.tolist()],
        dtype=np.int64,
    )
    if len(set(ids.tolist())) != len(ids):
        raise ValueError(f"Object Id values are not unique within region {region_name!r}")
    return ids


def _whole_number(value: object, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, (int, float, str, np.integer, np.floating)):
        raise ValueError(f"{label} is not a whole number: {value!r}")
    try:
        number = float(value)
    except ValueError as exc:
        raise ValueError(f"{label} is not a whole number: {value!r}") from exc
    if not np.isfinite(number):
        raise ValueError(f"{label} is not a whole number: {value!r}")
    rounded = round(number)
    if abs(number - rounded) > _WHOLE_NUMBER_TOLERANCE:
        raise ValueError(f"{label} is not a whole number: {value!r}")
    return int(rounded)


def _region_json(region_name: str, region: Mapping[str, object]) -> str:
    raw = region.get("json")
    if not isinstance(raw, str) or not raw:
        raise ValueError(f"Region {region_name!r} has no json segmentation path")
    return raw


def _region_image(regions: Mapping[str, object], region_name: str, region: Mapping[str, object]) -> str:
    avivator = regions.get("avivator")
    if not isinstance(avivator, dict) or not isinstance(avivator.get("base_url"), str):
        raise ValueError("regions.avivator.base_url is required to locate images")
    viv_image = region.get("viv_image")
    if not isinstance(viv_image, dict) or not isinstance(viv_image.get("file"), str) or not viv_image["file"]:
        raise ValueError(f"Region {region_name!r} has no viv_image.file")
    return str(Path(avivator["base_url"]) / viv_image["file"])


def _check_region_name(region_name: str) -> None:
    if not region_name or region_name in {".", ".."} or region_name != Path(region_name).name:
        raise ValueError(f"Region name {region_name!r} cannot be used as a store directory")


def _group_suffixes(group_name: str, columns: Sequence[str]) -> list[str] | None:
    prefix = group_name + "_"
    suffixes: list[str] = []
    for column in columns:
        if not column.startswith(prefix) or column[len(prefix) :] == "":
            return None
        suffixes.append(column[len(prefix) :])
    return suffixes


def _claim(consumed: set[str], columns: Sequence[str], group_name: str) -> None:
    overlap = [column for column in columns if column in consumed]
    if overlap:
        raise ValueError(
            f"Column group {group_name!r} uses columns already assigned: {', '.join(overlap)}"
        )
    consumed.update(columns)


def _columns_match(frame: pd.DataFrame, left: Sequence[str], right: Sequence[str]) -> bool:
    if len(left) != len(right):
        return False
    for source, target in zip(left, right):
        if source == target:
            continue
        source_values = np.asarray(pd.to_numeric(frame[source], errors="coerce"), dtype=np.float64)
        target_values = np.asarray(pd.to_numeric(frame[target], errors="coerce"), dtype=np.float64)
        if source_values.shape != target_values.shape or not np.allclose(
            source_values, target_values, equal_nan=True, rtol=0, atol=_WHOLE_NUMBER_TOLERANCE
        ):
            return False
    return True


def _columns_except_markers(frame: pd.DataFrame, layout: ColumnLayout) -> list[str]:
    marker_columns = {column for columns in layout.layers.values() for column in columns}
    return [column for column in frame.columns if column not in marker_columns]


def _groups_among(
    column_groups: Sequence[Mapping[str, object]],
    kept: set[str],
) -> list[tuple[str, list[str]]]:
    groups: list[tuple[str, list[str]]] = []
    for group in column_groups:
        name = str(group["name"])
        raw_columns = group.get("columns")
        if isinstance(raw_columns, str) or not isinstance(raw_columns, Sequence):
            continue
        columns = [str(column) for column in raw_columns]
        if columns and all(column in kept for column in columns):
            groups.append((name, columns))
    return groups


def _attach_regions(project: MDVProject, datasource: str, regions: Mapping[str, object]) -> None:
    import copy

    metadata = project.get_datasource_metadata(datasource)
    attached = copy.deepcopy(regions)
    all_regions = attached.get("all_regions")
    if isinstance(all_regions, dict):
        for name, region in all_regions.items():
            if isinstance(region, dict):
                region["spatial"] = {"coordinate_system": "global", "file": f"{name}.zarr"}
    metadata["regions"] = attached
    project.set_datasource_metadata(metadata)


def _symlink_project_dir(source: Path, output: Path, name: str) -> None:
    src = source / name
    if not src.is_dir():
        return
    os.symlink(src.resolve(), output / name, target_is_directory=True)


def _copy_views(source: MDVProject, project: MDVProject, mapping: Mapping[str, str]) -> None:
    views = source.views
    if not isinstance(views, dict) or not views:
        return
    state = source.state if isinstance(source.state, dict) else {}
    initial = state.get("initial_view")
    ordered = [name for name in (state.get("all_views") or []) if name in views]
    seen: set[str] = set()
    for name in ordered:
        project.set_view(
            name,
            rewrite_expression_fields(views[name], mapping),
            make_default=name == initial,
        )
        seen.add(name)
    for name, view in views.items():
        if name in seen:
            continue
        project.set_view(name, rewrite_expression_fields(view, mapping), make_default=name == initial)
