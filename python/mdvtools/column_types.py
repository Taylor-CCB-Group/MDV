"""The column-type contract (ADR-0024), as the Python package reads it.

A column's type has three parts: its stored type (a zarr v3 data type), its encoding (plain,
categorical or multi-valued) and its kind (numeric, categorical, unique or multi-valued). A
column's ``datatype`` is the legacy field.

The table of legacy datatypes is defined in ``src/datastore/columnTypes.ts``. ``column_types.json``
is a copy generated from it by the front end's tests, so both sides read one definition. Do not
edit the JSON by hand.
"""

import json
import re
from importlib import resources
from typing import Any

import numpy as np

_CONTRACT: dict[str, Any] = json.loads(
    resources.files("mdvtools").joinpath("column_types.json").read_text()
)

#: legacy datatype -> {storedType, encoding, kind, rowWidth}
LEGACY_DATATYPES: dict[str, dict[str, str]] = _CONTRACT["datatypes"]

_NUMPY_BY_STORED_TYPE: dict[str, type[np.generic]] = {
    "float32": np.float32,
    "int32": np.int32,
    "uint16": np.uint16,
    "uint8": np.uint8,
}

#: The numpy dtype each legacy datatype is stored as in ``datafile.h5``. ``unique`` is absent:
#: it is fixed-width bytes whose width depends on the values.
NUMPY_STORAGE_DTYPES: dict[str, type[np.generic]] = {
    datatype: _NUMPY_BY_STORED_TYPE[entry["storedType"]]
    for datatype, entry in LEGACY_DATATYPES.items()
    if entry["storedType"] in _NUMPY_BY_STORED_TYPE
}

_INTEGER = re.compile(r"u?int(8|16|32|64)")
_FLOAT = re.compile(r"float(16|32|64)")

#: numpy integer dtypes whose values all fit int32, so they can be stored exactly as `int32`.
#: Only non-nullable ones: int32 storage has no missing value, so pandas nullable integers
#: (`Int32`) and polars integers (any column may hold nulls) go to float32-backed `integer`.
_EXACT_INT32_NUMPY = frozenset({"int8", "int16", "int32", "uint8", "uint16"})

#: integer dtypes with values a float32 cannot hold exactly (above 2**24). Compared in lower case,
#: so this covers numpy, pandas nullable and polars names alike.
_WIDE_INTEGER_DTYPES = frozenset({"int32", "int64", "uint32", "uint64"})
_FLOAT32_EXACT_INTEGER_LIMIT = 2**24


def _is_polars(dtype: Any) -> bool:
    return type(dtype).__module__.startswith("polars") or getattr(dtype, "__module__", "").startswith(
        "polars"
    )


def _dtype_name(dtype: Any) -> str:
    """The dtype's name, case preserved: pandas nullable `Int32` differs from numpy `int32`."""
    # polars dtypes (classes or instances) print as e.g. "Int32" or "Datetime(time_unit='us', ...)";
    # numpy and pandas dtypes have a `name` such as "int32", "Int32" or "datetime64[ns, UTC]"
    if _is_polars(dtype):
        return str(dtype)
    return str(getattr(dtype, "name", dtype))


def infer_datatype(dtype: Any) -> str:
    """The legacy datatype a column of this numpy, pandas or polars dtype is stored as.

    Non-nullable numpy integers that fit int32 become ``int32``, stored exactly. Other integers
    (wider, pandas nullable, or polars) become ``integer``, which is float32: see
    :func:`integer_precision_warning`. Floats become ``double``. Dates become ``double`` (they are
    stored as days since the epoch). Everything else becomes ``text``, which may later turn into
    ``text16`` or ``unique`` depending on how many distinct values it has.
    """
    name = _dtype_name(dtype)
    lower = name.lower()
    if lower.startswith(("datetime", "date")):
        return "double"
    if _INTEGER.fullmatch(lower):
        if not _is_polars(dtype) and name in _EXACT_INT32_NUMPY:
            return "int32"
        return "integer"
    if _FLOAT.fullmatch(lower):
        return "double"
    return "text"


def integer_precision_warning(original_dtype: str, values: Any) -> str | None:
    """A ``storage_warnings`` entry if storing these values as float32 ``integer`` loses precision.

    Only integer dtypes wider than float32's exact range are checked. ``values`` may contain NaN.
    """
    if str(original_dtype).lower() not in _WIDE_INTEGER_DTYPES:
        return None
    numeric = np.asarray(values, dtype=np.float64)
    finite = np.abs(numeric[np.isfinite(numeric)])
    out_of_range = finite[finite > _FLOAT32_EXACT_INTEGER_LIMIT]
    if len(out_of_range) == 0:
        return None
    return (
        f"{original_dtype} -> float32 (integer): "
        f"{len(out_of_range)} value(s) exceed exact-integer range "
        f"(max abs={float(finite.max()):.0f}); precision loss possible."
    )
