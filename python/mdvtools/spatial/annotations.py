from __future__ import annotations

import logging
import os
import shutil
import warnings
from dataclasses import dataclass
from datetime import datetime
from typing import TYPE_CHECKING, Any

import h5py
import numpy
import pandas

from mdvtools.mdvproject import MDVProject

if TYPE_CHECKING:
    from anndata import AnnData


@dataclass
class SpatialAnnotationPatchReport:
    project_dir: str
    datasource: str
    spatialdata_path: str | None
    table_name: str
    region: str | list[str]
    region_key: str
    instance_key: str
    instance_count: int
    key_action: str
    preservation_column: str | None
    annotation_key: str
    annotation_columns: list[str]
    missing_annotation_count: int
    applied: bool
    backup_dir: str | None = None

    def to_text(self) -> str:
        lines = [
            "Spatial annotation patch report",
            f"- Project: {self.project_dir}",
            f"- Datasource: {self.datasource}",
            f"- SpatialData: {self.spatialdata_path or '<provided table>'}",
            f"- Table: {self.table_name}",
            f"- Region key: {self.region_key}",
            f"- Instance key: {self.instance_key}",
            f"- Rows: {self.instance_count}",
            f"- Key action: {self.key_action}",
            f"- Annotation key: {self.annotation_key}",
            f"- Annotation columns: {', '.join(self.annotation_columns) if self.annotation_columns else '<none>'}",
            f"- Missing annotations filled: {self.missing_annotation_count}",
            f"- Applied: {self.applied}",
        ]
        if self.preservation_column is not None:
            lines.append(f"- Preserved previous key values in: {self.preservation_column}")
        if self.backup_dir is not None:
            lines.append(f"- Backup: {self.backup_dir}")
        return "\n".join(lines)


def patch_spatial_annotations(
    project_dir: str,
    annotation_csv: str,
    *,
    datasource: str = "cells",
    spatialdata_path: str | None = None,
    table_name: str | None = None,
    annotation_key: str | None = None,
    separator: str = ",",
    missing_value: str = "ND",
    apply_changes: bool = False,
    overwrite_preserved: bool = False,
) -> SpatialAnnotationPatchReport:
    """Patch an existing MDV project using the SpatialData table instance key."""
    adata, resolved_spatialdata_path, resolved_table_name = _load_spatial_table(
        project_dir,
        spatialdata_path=spatialdata_path,
        table_name=table_name,
    )
    return patch_spatial_annotations_from_table(
        project_dir,
        annotation_csv,
        adata,
        table_name=resolved_table_name,
        spatialdata_path=resolved_spatialdata_path,
        datasource=datasource,
        annotation_key=annotation_key,
        separator=separator,
        missing_value=missing_value,
        apply_changes=apply_changes,
        overwrite_preserved=overwrite_preserved,
    )


def patch_spatial_annotations_from_table(
    project_dir: str,
    annotation_data: str | pandas.DataFrame,
    adata: "AnnData",
    *,
    table_name: str,
    spatialdata_path: str | None = None,
    datasource: str = "cells",
    annotation_key: str | None = None,
    separator: str = ",",
    missing_value: str = "ND",
    apply_changes: bool = False,
    overwrite_preserved: bool = False,
) -> SpatialAnnotationPatchReport:
    """Patch an existing MDV project from an already-loaded SpatialData table."""
    from spatialdata.models import get_table_keys

    _validate_project_dir(project_dir)

    try:
        region, region_key, instance_key = get_table_keys(adata)
    except ValueError as error:
        raise ValueError(
            "SpatialData table does not define table keys in uns['spatialdata_attrs']."
        ) from error

    if instance_key not in adata.obs.columns:
        raise ValueError(
            f"SpatialData instance key '{instance_key}' is not present in table obs columns."
        )

    source_ids = _coerce_key_values(adata.obs[instance_key], instance_key)
    _validate_no_duplicates(source_ids, f"SpatialData instance key '{instance_key}'")

    project = MDVProject(project_dir)
    ds_metadata = project.get_datasource_metadata(datasource)
    project_size = ds_metadata.get("size")
    if project_size != len(source_ids):
        raise ValueError(
            f"Datasource '{datasource}' has {project_size} rows but SpatialData table "
            f"'{table_name}' has {len(source_ids)} rows. V1 requires matching row order."
        )

    annotation_df = _load_annotation_dataframe(annotation_data, separator)
    annotation_key = annotation_key or instance_key
    annotation_df, annotation_columns = _normalize_annotation_dataframe(
        annotation_df,
        annotation_key=annotation_key,
        instance_key=instance_key,
    )
    missing_annotation_count = _validate_annotation_keys(
        annotation_df,
        annotation_key=instance_key,
        source_ids=source_ids,
    )

    current_fields = {column["field"] for column in ds_metadata["columns"]}
    current_key_values: list[Any] | None = None
    current_key_metadata = None
    key_action = "add"
    preservation_column: str | None = None

    if instance_key in current_fields:
        current_key_metadata = next(
            column for column in ds_metadata["columns"] if column["field"] == instance_key
        )
        current_key_values = _read_mdv_column(project_dir, ds_metadata, datasource, instance_key)
        current_key_strings = [str(value) for value in current_key_values]
        values_match = current_key_strings == source_ids
        if values_match and current_key_metadata.get("datatype") == "unique":
            key_action = "unchanged"
        elif values_match:
            key_action = "rewrite datatype"
        else:
            key_action = "repair"
            preservation_column = _preservation_column(instance_key)
            if preservation_column in current_fields and not overwrite_preserved:
                raise ValueError(
                    f"Cannot preserve previous '{instance_key}' values because "
                    f"'{preservation_column}' already exists in datasource '{datasource}'. "
                    "Pass overwrite_preserved=True to replace it."
                )

    if apply_changes:
        backup_dir = _backup_project_files(project_dir)
        if preservation_column is not None:
            assert current_key_values is not None
            project.set_column(
                datasource,
                {
                    "name": preservation_column,
                    "field": preservation_column,
                    "datatype": "unique",
                },
                [str(value) for value in current_key_values],
            )
        if key_action != "unchanged":
            project.set_column(
                datasource,
                {
                    "name": instance_key,
                    "field": instance_key,
                    "datatype": "unique",
                },
                source_ids,
            )
        _write_annotation_columns(
            project,
            datasource=datasource,
            annotation_df=annotation_df,
            annotation_columns=annotation_columns,
            annotation_key=instance_key,
            source_ids=source_ids,
            missing_value=missing_value,
        )
    else:
        backup_dir = None

    return SpatialAnnotationPatchReport(
        project_dir=project_dir,
        datasource=datasource,
        spatialdata_path=spatialdata_path,
        table_name=table_name,
        region=region,
        region_key=region_key,
        instance_key=instance_key,
        instance_count=len(source_ids),
        key_action=key_action,
        preservation_column=preservation_column,
        annotation_key=annotation_key,
        annotation_columns=annotation_columns,
        missing_annotation_count=missing_annotation_count,
        applied=apply_changes,
        backup_dir=backup_dir,
    )


def _validate_project_dir(project_dir: str) -> None:
    required_files = ("datafile.h5", "datasources.json")
    missing = [
        filename for filename in required_files
        if not os.path.exists(os.path.join(project_dir, filename))
    ]
    if missing:
        raise FileNotFoundError(
            f"'{project_dir}' does not appear to be an MDV project "
            f"(missing {', '.join(missing)})."
        )


def _load_spatial_table(
    project_dir: str,
    *,
    spatialdata_path: str | None,
    table_name: str | None,
) -> tuple["AnnData", str, str]:
    from anndata import AnnData

    resolved_path = spatialdata_path or _detect_single_spatialdata_path(project_dir)
    sdata = _read_spatialdata_zarr(resolved_path)

    table_items = list(sdata.tables.items())
    if table_name is None:
        if len(table_items) != 1:
            table_names = ", ".join(name for name, _adata in table_items) or "<none>"
            raise ValueError(
                "V1 requires exactly one SpatialData table when --table-name is not supplied. "
                f"Found {len(table_items)} table(s): {table_names}."
            )
        resolved_table_name, adata = table_items[0]
    else:
        if table_name not in sdata.tables:
            table_names = ", ".join(name for name, _adata in table_items) or "<none>"
            raise ValueError(
                f"SpatialData table '{table_name}' was not found. Available tables: {table_names}."
            )
        resolved_table_name = table_name
        adata = sdata.tables[table_name]

    if not isinstance(adata, AnnData):
        raise ValueError(f"SpatialData table '{resolved_table_name}' is not an AnnData table.")
    return adata, resolved_path, resolved_table_name


def _detect_single_spatialdata_path(project_dir: str) -> str:
    spatial_dir = os.path.join(project_dir, "spatial")
    if not os.path.isdir(spatial_dir):
        raise FileNotFoundError(
            f"Cannot auto-detect SpatialData because '{spatial_dir}' does not exist."
        )

    candidates: list[str] = []
    for path in [spatial_dir] + [
        os.path.join(spatial_dir, child)
        for child in sorted(os.listdir(spatial_dir))
        if not child.startswith(".")
    ]:
        if not os.path.isdir(path):
            continue
        try:
            _read_spatialdata_zarr(path)
        except Exception:
            continue
        candidates.append(path)

    if len(candidates) != 1:
        candidate_text = ", ".join(candidates) if candidates else "<none>"
        raise ValueError(
            "V1 requires exactly one SpatialData store under the project spatial folder. "
            f"Found {len(candidates)} candidate(s): {candidate_text}."
        )
    return candidates[0]


def _read_spatialdata_zarr(path: str):
    import spatialdata as sd

    ome_logger = logging.getLogger("ome_zarr.reader")
    original_level = ome_logger.level
    ome_logger.setLevel(logging.ERROR)
    try:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=UserWarning)
            return sd.read_zarr(path)
    finally:
        ome_logger.setLevel(original_level)


def _coerce_key_values(values: Any, key_name: str) -> list[str]:
    if not isinstance(values, pandas.Series):
        values = pandas.Series(values)
    if values.isna().any():
        raise ValueError(f"Key column '{key_name}' contains missing values.")
    return [str(value) for value in values.tolist()]


def _validate_no_duplicates(values: list[str], label: str) -> None:
    seen: set[str] = set()
    duplicates: list[str] = []
    for value in values:
        if value in seen and value not in duplicates:
            duplicates.append(value)
        seen.add(value)
    if duplicates:
        examples = ", ".join(duplicates[:5])
        raise ValueError(
            f"{label} contains {len(duplicates)} duplicate value(s), "
            f"for example: {examples}."
        )


def _load_annotation_dataframe(
    annotation_data: str | pandas.DataFrame,
    separator: str,
) -> pandas.DataFrame:
    if isinstance(annotation_data, str):
        return pandas.read_csv(annotation_data, sep=separator)
    return annotation_data.copy()


def _normalize_annotation_dataframe(
    annotation_df: pandas.DataFrame,
    *,
    annotation_key: str,
    instance_key: str,
) -> tuple[pandas.DataFrame, list[str]]:
    if annotation_key not in annotation_df.columns:
        raise ValueError(
            f"Annotation key '{annotation_key}' was not found in the annotation data."
        )
    if annotation_key != instance_key and instance_key in annotation_df.columns:
        raise ValueError(
            f"Cannot rename annotation key '{annotation_key}' to '{instance_key}' "
            f"because the annotation data already contains a '{instance_key}' column."
        )
    if annotation_key != instance_key:
        annotation_df = annotation_df.rename(columns={annotation_key: instance_key})
    annotation_columns = [
        column for column in annotation_df.columns
        if column != instance_key
    ]
    return annotation_df, annotation_columns


def _validate_annotation_keys(
    annotation_df: pandas.DataFrame,
    *,
    annotation_key: str,
    source_ids: list[str],
) -> int:
    annotation_ids = _coerce_key_values(annotation_df[annotation_key], annotation_key)
    _validate_no_duplicates(annotation_ids, f"Annotation key '{annotation_key}'")
    annotation_df[annotation_key] = annotation_ids

    source_id_set = set(source_ids)
    annotation_id_set = set(annotation_ids)
    extra_ids = annotation_id_set - source_id_set
    if extra_ids:
        examples = ", ".join(sorted(extra_ids)[:5])
        raise ValueError(
            f"Annotation data contains {len(extra_ids)} id(s) not present in the "
            f"SpatialData instance key '{annotation_key}', for example: {examples}."
        )
    return len(source_id_set - annotation_id_set)


def _preservation_column(instance_key: str) -> str:
    if instance_key == "cell_id":
        return "cell_index"
    return f"{instance_key}_index"


def _backup_project_files(project_dir: str) -> str:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    backup_dir = os.path.join(project_dir, "annotation_patch_backups", timestamp)
    os.makedirs(backup_dir, exist_ok=False)
    for filename in ("datafile.h5", "datasources.json"):
        shutil.copy2(os.path.join(project_dir, filename), os.path.join(backup_dir, filename))
    return backup_dir


def _read_mdv_column(
    project_dir: str,
    datasource_metadata: dict[str, Any],
    datasource: str,
    column: str,
) -> list[Any]:
    column_metadata = next(
        metadata for metadata in datasource_metadata["columns"]
        if metadata["field"] == column
    )
    with h5py.File(os.path.join(project_dir, "datafile.h5"), "r") as h5:
        group = h5[datasource]
        if not isinstance(group, h5py.Group):
            raise AttributeError(f"'{datasource}' is not a group in datafile.h5")
        raw_data = numpy.array(group[column])

    datatype = column_metadata["datatype"]
    if datatype in {"text", "text16"}:
        values = column_metadata["values"]
        return [values[index] for index in raw_data]
    if datatype == "multitext":
        string_length = column_metadata["stringLength"]
        chunks = numpy.split(raw_data, int(raw_data.shape[0] / string_length))
        values = column_metadata["values"]
        return [
            ",".join(values[index] for index in chunk if index != 65535)
            for chunk in chunks
        ]
    if datatype == "unique":
        return [
            value.decode("utf-8") if isinstance(value, bytes) else str(value)
            for value in raw_data
        ]
    return list(raw_data)


def _write_annotation_columns(
    project: MDVProject,
    *,
    datasource: str,
    annotation_df: pandas.DataFrame,
    annotation_columns: list[str],
    annotation_key: str,
    source_ids: list[str],
    missing_value: str,
) -> None:
    annotation_by_id = annotation_df.set_index(annotation_key)
    for column in annotation_columns:
        values_by_id = annotation_by_id[column].to_dict()
        values = [
            values_by_id.get(source_id, missing_value)
            for source_id in source_ids
        ]
        values = [
            missing_value if pandas.isna(value) else value
            for value in values
        ]
        project.set_column(
            datasource,
            {
                "name": column,
                "field": column,
                "datatype": "text",
            },
            values,
        )
