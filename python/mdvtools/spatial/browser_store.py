"""Browser-oriented rewrites of the SpatialData store copied into a project.

`convert-spatial` copies each source store into ``<project>/spatial``. These options
have spatialdata-js-util rewrite that copy so SpatialData.js can read it cheaply over
HTTP: images recompressed as HTJ2K or JPEG 2000 (with a pyramid for single-level
images), table matrices in CSC, points Morton-sorted, and an experimental points-zarr
copy beside each points element.

Presets bundle a recommended set and are expected to change as the formats and the
frontend evolve; the per-item options override whatever the preset chose. The
rewrites only ever touch the project's copy, so they cannot be combined with linking.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable

from .conversion_options import DEFAULT_BROWSER_PRESET, IMAGE_PRESET_CHOICES

if TYPE_CHECKING:
    from spatialdata import SpatialData


@dataclass(frozen=True)
class BrowserStoreOptions:
    # "htj2k" or "jp2k"; None leaves images as spatialdata writes them.
    image_codec: str | None = None
    image_preset: str = "balanced"
    csc_tables: bool = False
    morton_points: bool = False
    points_zarr: bool = False

    def __post_init__(self) -> None:
        if self.image_codec is not None and self.image_codec not in _CODEC_NAMES:
            raise ValueError(
                f"Unknown image codec '{self.image_codec}'. Expected one of: {', '.join(_CODEC_NAMES)}"
            )
        if self.image_preset not in IMAGE_PRESET_CHOICES:
            raise ValueError(
                f"Unknown image preset '{self.image_preset}'. "
                f"Expected one of: {', '.join(IMAGE_PRESET_CHOICES)}"
            )

    @property
    def enabled(self) -> bool:
        return (
            self.image_codec is not None
            or self.csc_tables
            or self.morton_points
            or self.points_zarr
        )

    def describe(self) -> str:
        parts = []
        if self.image_codec is not None:
            parts.append(f"{self.image_codec} images ({self.image_preset})")
        if self.csc_tables:
            parts.append("CSC tables")
        if self.morton_points:
            parts.append("Morton-sorted points")
        if self.points_zarr:
            parts.append("points-zarr copies")
        return ", ".join(parts) if parts else "none"


_CODEC_NAMES = {
    "htj2k": "experimental.openjph_htj2k",
    "jp2k": "imagecodecs_jpeg2k",
}

_STANDARD = BrowserStoreOptions(
    image_codec="htj2k",
    image_preset="balanced",
    csc_tables=True,
    morton_points=True,
)

BROWSER_PRESETS: dict[str, BrowserStoreOptions] = {
    "none": BrowserStoreOptions(),
    "standard": _STANDARD,
    # points-zarr is a draft spec that published SpatialData.js does not read yet.
    "experimental": replace(_STANDARD, points_zarr=True),
}


def resolve_browser_store_options(
    preset: str = DEFAULT_BROWSER_PRESET,
    *,
    image_codec: str | None = None,
    image_preset: str | None = None,
    csc_tables: bool | None = None,
    morton_points: bool | None = None,
    points_zarr: bool | None = None,
) -> BrowserStoreOptions:
    """Start from *preset* and apply every override that is not None.

    *image_codec* ``"none"`` turns image recompression off even if the preset has it on.
    """
    if preset not in BROWSER_PRESETS:
        raise ValueError(
            f"Unknown browser preset '{preset}'. Expected one of: {', '.join(BROWSER_PRESETS)}"
        )
    overrides: dict[str, Any] = {}
    if image_codec is not None:
        overrides["image_codec"] = None if image_codec == "none" else image_codec
    if image_preset is not None:
        overrides["image_preset"] = image_preset
    if csc_tables is not None:
        overrides["csc_tables"] = csc_tables
    if morton_points is not None:
        overrides["morton_points"] = morton_points
    if points_zarr is not None:
        overrides["points_zarr"] = points_zarr
    return replace(BROWSER_PRESETS[preset], **overrides)


def write_spatialdata_store(
    sdata: SpatialData,
    dest: str | Path,
    options: BrowserStoreOptions,
    progress: Callable[[str], None] = print,
) -> None:
    """Write *sdata* to *dest*, then apply the rewrites *options* asks for in place."""
    dest = Path(dest)
    image_keys = _recompressible_image_keys(sdata) if options.image_codec is not None else []
    if options.image_codec is not None:
        for key in sdata.images:
            if key not in image_keys:
                progress(
                    f"Not recompressing images/{key}: browser codecs need integer pixels of 16 bits or fewer"
                )
    if options.image_codec is None or not image_keys:
        sdata.write(dest)
    else:
        from spatialdata_js_util import recompress_spatialdata

        progress(f"Writing '{dest.name}' with {options.image_codec} images ({options.image_preset})")
        # Writes the store, then recompresses the listed rasters; single-level images
        # get a pyramid first, existing pyramids are kept.
        recompress_spatialdata(
            sdata,
            dest,
            config={"images": {key: {} for key in image_keys}},
            codec=_CODEC_NAMES[options.image_codec],
            preset=_image_preset(options.image_preset),
            pyramid=True,
            manifest=False,
        )

    if options.csc_tables:
        from spatialdata_js_util import convert_store_tables_to_csc

        progress(f"Converting tables in '{dest.name}' to CSC")
        convert_store_tables_to_csc(dest, manifest=False)

    if options.points_zarr or options.morton_points:
        import pyarrow.dataset as ds
        from spatialdata_js_util import MORTON_CODE_2D_COLUMN, add_points_zarr_to_store
        from spatialdata_js_util.runners import run_morton_points_from_zarr
        from spatialdata_js_util.store import (
            list_points_keys,
            points_parquet_path,
            read_points_element_attrs,
        )

        for points_key in list_points_keys(dest):
            if options.points_zarr:
                if read_points_element_attrs(dest, points_key)["feature_key"] is None:
                    progress(f"No points-zarr copy of points/{points_key}: it declares no feature_key")
                else:
                    progress(f"Adding points-zarr copy of points/{points_key}")
                    add_points_zarr_to_store(None, dest, points_key=points_key, progress=progress)
            if options.morton_points:
                schema = ds.dataset(points_parquet_path(dest, points_key), format="parquet").schema
                if MORTON_CODE_2D_COLUMN in schema.names:
                    # Already Morton-indexed, possibly in a deliberate order such as
                    # feature-then-Morton that a plain re-sort would undo.
                    progress(f"Keeping points/{points_key}: already Morton-indexed")
                else:
                    progress(f"Morton-sorting points/{points_key}")
                    run_morton_points_from_zarr(dest, points_key=points_key)


def _recompressible_image_keys(sdata: SpatialData) -> list[str]:
    from spatialdata_js_util.images import SUPPORTED_BROWSER_JP2K_DTYPES
    from xarray import DataTree

    keys = []
    for key, element in sdata.images.items():
        # A multiscale image is a DataTree of scale nodes, each holding an "image" variable.
        raster = next(iter(element.values()))["image"] if isinstance(element, DataTree) else element
        if raster.dtype in SUPPORTED_BROWSER_JP2K_DTYPES:
            keys.append(key)
    return keys


def _image_preset(preset: str):
    # Narrows to sd-js-util's ImagePreset literal without a cast; __post_init__
    # has already rejected anything else.
    if preset == "lossless":
        return "lossless"
    if preset == "balanced":
        return "balanced"
    if preset == "small":
        return "small"
    raise ValueError(f"Unknown image preset '{preset}'")
