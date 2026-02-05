"""
Helper functions for merging annotation anndata objects into SpatialData tables.

This module provides utilities for merging additional anndata objects (e.g., containing
UMAP annotations) into tables extracted from SpatialData objects during conversion.

Current scope:
- Only supports merging obs columns onto existing instances (mapped via instance_key)
- Merges obsm (dimensionality reductions) from annotation data
- Skips X layer from annotation (keeps original table's X)

Future considerations (not currently supported):
- Adding layers: Would require handling multiple expression matrices per table
- Multi-modal data (mu-data like mappings): Would require either:
  * Multiple tables with redundant 'obs' (one per modality)
  * Single table where some 'vars' refer to 'genes' and others to 'rna', 'atac', etc.
  This would need significant changes to how tables are structured and how MDV handles
  gene/feature identifiers.
"""

from typing import TYPE_CHECKING, Optional, Tuple
import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from anndata import AnnData

def _infer_merge_key(
    original_adata: "AnnData",
    annotation_adata: "AnnData",
    user_specified_key: Optional[str] = None
) -> Tuple[str, dict]:
    """
    Infer the merge key column from both anndata objects using get_table_keys.
    
    This function requires that we can identify an instance_key in the original table
    and map annotation data to it. We error out early if this isn't possible to avoid
    writing bad data to disk.
    
    Args:
        original_adata: The original anndata from SpatialData table
        annotation_adata: The annotation anndata to merge in
        user_specified_key: Optional user-specified merge key column name
        
    Returns:
        Tuple of (merge_key_column_name, diagnostic_stats_dict)
        
    Raises:
        ValueError: If no suitable merge key can be found or instance_key cannot be determined
    """
    from spatialdata.models import get_table_keys
    
    diagnostic_stats = {
        "original_instance_key": None,
        "annotation_instance_key": None,
        "merge_key_found": False,
        "merge_key_source": None,
    }
    
    # First, try to get instance_key from original table - this is required
    try:
        _r_orig, _obs_region_col_orig, instance_key_orig = get_table_keys(original_adata)
        diagnostic_stats["original_instance_key"] = instance_key_orig
    except Exception as e:
        raise ValueError(
            f"Cannot merge annotation data: failed to get table keys from original anndata: {e}. "
            f"Merge requires that the original table has a valid instance_key."
        ) from e
    
    if instance_key_orig is None:
        raise ValueError(
            "Cannot merge annotation data: original table has no instance_key. "
            "Merge requires mapping annotation data to existing instances via instance_key."
        )
    
    if instance_key_orig not in original_adata.obs.columns:
        raise ValueError(
            f"Cannot merge annotation data: instance_key '{instance_key_orig}' from get_table_keys "
            f"not found in original table's obs columns. Available columns: {list(original_adata.obs.columns)}"
        )
    
    # Try to get instance_key from annotation data (preferred but not strictly required)
    try:
        _r_ann, _obs_region_col_ann, instance_key_ann = get_table_keys(annotation_adata)
        diagnostic_stats["annotation_instance_key"] = instance_key_ann
    except Exception as e:
        # This is a warning, not an error - we can still merge if we find a matching column
        print(f"WARNING: Could not get table keys from annotation anndata: {e}")
        instance_key_ann = None
    
    # If user specified a key, validate it exists and can map
    if user_specified_key:
        if user_specified_key not in original_adata.obs.columns:
            raise ValueError(
                f"Cannot merge: user-specified merge key '{user_specified_key}' not found in original table. "
                f"Available columns: {list(original_adata.obs.columns)}"
            )
        if user_specified_key not in annotation_adata.obs.columns:
            raise ValueError(
                f"Cannot merge: user-specified merge key '{user_specified_key}' not found in annotation table. "
                f"Available columns: {list(annotation_adata.obs.columns)}"
            )
        diagnostic_stats["merge_key_found"] = True
        diagnostic_stats["merge_key_source"] = "user_specified"
        return user_specified_key, diagnostic_stats
    
    # Prefer instance_key from original table if it exists in annotation
    if instance_key_orig in annotation_adata.obs.columns:
        diagnostic_stats["merge_key_found"] = True
        diagnostic_stats["merge_key_source"] = "instance_key_match"
        return instance_key_orig, diagnostic_stats
    
    # If annotation has instance_key, try that if it exists in original
    if instance_key_ann and instance_key_ann in original_adata.obs.columns:
        diagnostic_stats["merge_key_found"] = True
        diagnostic_stats["merge_key_source"] = "annotation_instance_key"
        return instance_key_ann, diagnostic_stats
    
    # Fallback: try common column names that might map to instance_key
    common_columns = set(original_adata.obs.columns) & set(annotation_adata.obs.columns)
    
    # Priority order for fallback (common names that might be instance identifiers)
    fallback_names = ["EntityID", "entity_id", "cell_id", "CellID"]
    
    for fallback_name in fallback_names:
        if fallback_name in common_columns:
            print(f"INFO: Using fallback merge key '{fallback_name}' (instance_key not directly mappable).")
            diagnostic_stats["merge_key_found"] = True
            diagnostic_stats["merge_key_source"] = f"fallback_{fallback_name}"
            return fallback_name, diagnostic_stats
    
    # No suitable key found - error out
    raise ValueError(
        f"Cannot merge annotation data: no suitable merge key found. "
        f"Original table instance_key: '{instance_key_orig}' (not in annotation). "
        f"Annotation instance_key: '{instance_key_ann}' (not in original). "
        f"Common columns: {sorted(common_columns)}. "
        f"Merge requires that annotation data can map to at least one instance_key in the original table."
    )


def _convert_merge_key_dtype(
    original_adata: "AnnData",
    annotation_adata: "AnnData",
    merge_key: str
) -> Tuple["AnnData", "AnnData", list[str]]:
    """
    Convert merge key columns to compatible dtypes (prefer string/object).
    
    Args:
        original_adata: Original anndata object
        annotation_adata: Annotation anndata object
        merge_key: Name of the merge key column
        
    Returns:
        Tuple of (original_adata_copy, annotation_adata_copy, warnings)
    """
    warnings = []
    
    # Make copies to avoid mutating originals
    original_adata = original_adata.copy()
    annotation_adata = annotation_adata.copy()
    
    orig_dtype = original_adata.obs[merge_key].dtype
    ann_dtype = annotation_adata.obs[merge_key].dtype
    
    # If dtypes already match, no conversion needed
    if orig_dtype == ann_dtype:
        return original_adata, annotation_adata, warnings
    
    # Convert both to string/object for compatibility
    if orig_dtype != 'object':
        warnings.append(
            f"INFO: Converting merge key '{merge_key}' in original anndata from {orig_dtype} to object/string."
        )
        original_adata.obs[merge_key] = original_adata.obs[merge_key].astype(str)
    
    if ann_dtype != 'object':
        warnings.append(
            f"INFO: Converting merge key '{merge_key}' in annotation anndata from {ann_dtype} to object/string."
        )
        annotation_adata.obs[merge_key] = annotation_adata.obs[merge_key].astype(str)
    
    return original_adata, annotation_adata, warnings


def _compute_merge_diagnostics(
    original_adata: "AnnData",
    annotation_adata: "AnnData",
    merge_key: str
) -> dict:
    """
    Compute diagnostic statistics about how well the annotation data maps to original instances.
    
    Args:
        original_adata: Original anndata object
        annotation_adata: Annotation anndata object
        merge_key: Name of the merge key column
        
    Returns:
        Dictionary with diagnostic statistics
    """
    # Get unique values
    orig_keys = set(original_adata.obs[merge_key].dropna().unique())
    ann_keys = set(annotation_adata.obs[merge_key].dropna().unique())
    
    # Compute overlap statistics
    overlap = orig_keys & ann_keys
    orig_only = orig_keys - ann_keys
    ann_only = ann_keys - orig_keys
    
    overlap_pct = len(overlap) / len(orig_keys) * 100 if len(orig_keys) > 0 else 0
    coverage_pct = len(overlap) / len(ann_keys) * 100 if len(ann_keys) > 0 else 0
    
    # Check for duplicates
    orig_dups = original_adata.obs[merge_key].duplicated().sum()
    ann_dups = annotation_adata.obs[merge_key].duplicated().sum()
    
    diagnostics = {
        "original_unique_keys": len(orig_keys),
        "annotation_unique_keys": len(ann_keys),
        "overlapping_keys": len(overlap),
        "original_only_keys": len(orig_only),
        "annotation_only_keys": len(ann_only),
        "overlap_percentage": overlap_pct,
        "coverage_percentage": coverage_pct,
        "original_duplicates": int(orig_dups),
        "annotation_duplicates": int(ann_dups),
    }
    
    return diagnostics


def _merge_annotation_anndata(
    original_adata: "AnnData",
    annotation_adata: "AnnData",
    merge_key_column: Optional[str] = None,
    annotation_path: Optional[str] = None
) -> Tuple["AnnData", list[str], dict]:
    """
    Merge annotation anndata into original anndata.
    
    Current scope: Only merges obs columns onto existing instances (mapped via instance_key).
    Also merges obsm (dimensionality reductions) from annotation data.
    Skips X layer from annotation (keeps original table's X).
    
    Future considerations (not currently supported):
    - Adding layers: Would require handling multiple expression matrices per table
    - Multi-modal data: Would require restructuring tables to handle vars that refer to
      different modalities (genes, rna, atac, etc.)
    
    Args:
        original_adata: The original anndata from SpatialData table
        annotation_adata: The annotation anndata to merge in
        merge_key_column: Optional user-specified merge key column name
        annotation_path: Optional path to annotation anndata file (for provenance)
        
    Returns:
        Tuple of (merged_anndata, list_of_warnings, provenance_dict)
        If merge is not possible (e.g., no instance_key), returns original_adata unchanged
        with provenance indicating failure.
        
    Note:
        This function does not raise exceptions - it returns the original anndata unchanged
        if merge is not possible. The calling code should check if at least one table was
        successfully merged across all tables in the SpatialData object.
    """
    import datetime
    
    all_warnings = []
    
    # Try to infer merge key - if this fails, we return original unchanged
    try:
        merge_key, key_diagnostics = _infer_merge_key(original_adata, annotation_adata, merge_key_column)
    except ValueError as e:
        # Merge not possible for this table - return original unchanged with failure provenance
        provenance = {
            "annotation_path": annotation_path,
            "merge_timestamp": datetime.datetime.now().isoformat(),
            "merge_key": None,
            "merge_key_source": None,
            "original_n_obs": original_adata.n_obs,
            "original_n_vars": original_adata.n_vars,
            "annotation_n_obs": annotation_adata.n_obs,
            "annotation_n_vars": annotation_adata.n_vars,
            "merged_obs_columns": [],
            "merged_obsm_keys": [],
            "merge_diagnostics": {},
            "merge_successful": False,
            "merge_failure_reason": str(e),
            "warnings": [f"Merge skipped: {str(e)}"],
        }
        original_adata.uns.setdefault("mdv", {})
        original_adata.uns["mdv"]["merge_provenance"] = provenance
        all_warnings.append(f"WARNING: Cannot merge annotation data into this table: {e}")
        return original_adata, all_warnings, provenance
    
    # Compute diagnostic statistics about fit quality
    merge_diagnostics = _compute_merge_diagnostics(original_adata, annotation_adata, merge_key)
    
    # Log diagnostic stats
    print(f"\n=== Annotation merge diagnostics ===")
    print(f"Merge key: '{merge_key}' (source: {key_diagnostics['merge_key_source']})")
    print(f"Original instances: {merge_diagnostics['original_unique_keys']}")
    print(f"Annotation instances: {merge_diagnostics['annotation_unique_keys']}")
    print(f"Overlapping instances: {merge_diagnostics['overlapping_keys']} ({merge_diagnostics['overlap_percentage']:.1f}% of original)")
    print(f"Original-only instances: {merge_diagnostics['original_only_keys']}")
    print(f"Annotation-only instances: {merge_diagnostics['annotation_only_keys']}")
    if merge_diagnostics['original_duplicates'] > 0:
        print(f"WARNING: {merge_diagnostics['original_duplicates']} duplicate keys in original")
    if merge_diagnostics['annotation_duplicates'] > 0:
        print(f"WARNING: {merge_diagnostics['annotation_duplicates']} duplicate keys in annotation")
    print(f"=====================================\n")
    
    # Warn if overlap is low
    if merge_diagnostics['overlap_percentage'] < 50:
        all_warnings.append(
            f"WARNING: Only {merge_diagnostics['overlap_percentage']:.1f}% of original instances have matches in annotation data. "
            f"Merge will result in many missing annotation values."
        )
    elif merge_diagnostics['overlap_percentage'] < 90:
        all_warnings.append(
            f"INFO: {merge_diagnostics['overlap_percentage']:.1f}% of original instances have matches in annotation data. "
            f"Some rows will have missing annotation values."
        )
    
    if merge_diagnostics['original_only_keys'] > 0:
        all_warnings.append(
            f"INFO: {merge_diagnostics['original_only_keys']} instances in original data have no match in annotation "
            f"(will keep original values for these rows)."
        )
    
    if merge_diagnostics['annotation_only_keys'] > 0:
        all_warnings.append(
            f"INFO: {merge_diagnostics['annotation_only_keys']} instances in annotation data have no match in original "
            f"(these will be skipped in merge)."
        )
    
    # Initialize provenance tracking
    provenance = {
        "annotation_path": annotation_path,
        "merge_timestamp": datetime.datetime.now().isoformat(),
        "merge_key": merge_key,
        "merge_key_source": key_diagnostics['merge_key_source'],
        "original_n_obs": original_adata.n_obs,
        "original_n_vars": original_adata.n_vars,
        "annotation_n_obs": annotation_adata.n_obs,
        "annotation_n_vars": annotation_adata.n_vars,
        "merged_obs_columns": [],
        "merged_obsm_keys": [],
        "merge_diagnostics": merge_diagnostics,
        "warnings": all_warnings.copy(),
    }
    
    # Make copies to avoid mutating originals
    merged_adata = original_adata.copy()
    annotation_adata_copy = annotation_adata.copy()
    
    # Convert merge key dtypes if needed
    merged_adata, annotation_adata_copy, dtype_warnings = _convert_merge_key_dtype(
        merged_adata, annotation_adata_copy, merge_key
    )
    all_warnings.extend(dtype_warnings)
    
    # Print warnings
    for warning in all_warnings:
        print(warning)
    
    # Merge obs columns using key-based merge
    # Strategy: prefer annotation values, keep original if annotation is missing/NaN
    # Set merge key as index for both
    merged_obs = merged_adata.obs.set_index(merge_key)
    ann_obs = annotation_adata_copy.obs.set_index(merge_key)
    
    # Merge with left join (keep all original rows)
    # For matching columns, prefer annotation values but keep original where annotation is missing
    merged_obs = merged_obs.join(ann_obs, how='left', rsuffix='_ann')
    
    # Handle columns that exist in both
    for col in annotation_adata_copy.obs.columns:
        if col == merge_key:
            continue
        ann_col = col
        if ann_col in merged_obs.columns and f"{ann_col}_ann" in merged_obs.columns:
            # Column exists in both - prefer annotation, fallback to original
            ann_values = merged_obs[f"{ann_col}_ann"]
            # Use annotation where available, otherwise original
            mask = pd.notna(ann_values)
            merged_obs.loc[mask, ann_col] = ann_values[mask]
            # Drop the _ann suffix column
            merged_obs = merged_obs.drop(columns=[f"{ann_col}_ann"])
            provenance["merged_obs_columns"].append({"name": col, "source": "merged"})
        elif f"{ann_col}_ann" in merged_obs.columns:
            # New column from annotation (was _ann suffix)
            merged_obs[ann_col] = merged_obs[f"{ann_col}_ann"]
            merged_obs = merged_obs.drop(columns=[f"{ann_col}_ann"])
            provenance["merged_obs_columns"].append({"name": col, "source": "new"})
    
    # Reset index to restore merge_key as column
    merged_adata.obs = merged_obs.reset_index()
    
    # Merge obsm (dimensionality reductions like UMAP)
    # Annotation takes precedence for matching keys
    for key in annotation_adata_copy.obsm.keys():
        if key in merged_adata.obsm:
            all_warnings.append(f"INFO: Replacing obsm['{key}'] with annotation version.")
        merged_adata.obsm[key] = annotation_adata_copy.obsm[key].copy()
        provenance["merged_obsm_keys"].append(key)
    
    # Note: We do NOT merge var or varm - current scope is obs-only merging
    # Future: If we want to support multi-modal data where vars refer to different
    # modalities (genes, rna, atac, etc.), we would need to restructure how tables
    # handle feature identifiers. This is out of scope for now.
    
    # Note: We intentionally skip X layer from annotation - only use original table's X
    # This is to avoid redundant expression matrices
    provenance["x_layer_skipped"] = True
    
    # Merge uns metadata (annotation takes precedence for matching keys)
    for key in annotation_adata_copy.uns.keys():
        if key in merged_adata.uns and key != "mdv":  # Preserve mdv metadata from original
            all_warnings.append(f"INFO: Replacing uns['{key}'] with annotation version.")
        if key != "mdv":  # Don't overwrite mdv metadata
            merged_adata.uns[key] = annotation_adata_copy.uns[key]
    
    # Store provenance in uns["mdv"]["merge_provenance"]
    merged_adata.uns.setdefault("mdv", {})
    provenance["merge_successful"] = True
    provenance["final_n_obs"] = merged_adata.n_obs
    provenance["final_n_vars"] = merged_adata.n_vars
    merged_adata.uns["mdv"]["merge_provenance"] = provenance
    
    return merged_adata, all_warnings, provenance

