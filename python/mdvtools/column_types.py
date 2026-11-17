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


def _dtype_name(dtype: Any) -> str:
    # polars dtypes (classes or instances) print as e.g. "Int32" or "Datetime(time_unit='us', ...)";
    # numpy and pandas dtypes have a `name` such as "int32", "Int32" or "datetime64[ns, UTC]"
    if type(dtype).__module__.startswith("polars") or getattr(dtype, "__module__", "").startswith(
        "polars"
    ):
        return str(dtype).lower()
    return str(getattr(dtype, "name", dtype)).lower()


def infer_datatype(dtype: Any) -> str:
    """The legacy datatype a column of this numpy, pandas or polars dtype is stored as.

    Integers become ``integer`` and floats ``double``, nullable or not. Dates become ``double``
    (they are stored as days since the epoch). Everything else becomes ``text``, which may later
    turn into ``text16`` or ``unique`` depending on how many distinct values it has.
    """
    name = _dtype_name(dtype)
    if name.startswith(("datetime", "date")):
        return "double"
    if _INTEGER.fullmatch(name):
        return "integer"
    if _FLOAT.fullmatch(name):
        return "double"
    return "text"
