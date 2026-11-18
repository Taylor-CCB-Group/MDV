# One column-type contract, built on zarr data types

**Status:** proposed. Addresses #480. Comes before ADR-0023's upgrade script and zarr loader, so
that neither adds yet another type mapping.

## Context

A column's `datatype` mixes how its values are stored with how charts treat them, and the mapping
from source dtypes is written several times over and disagrees with itself.

On the Python side:

- `datatype_mappings` ([mdvproject.py:149](../../python/mdvtools/mdvproject.py)) maps numpy
  `int32` to `double`, and `bool` to text. Nullable pandas dtypes such as `Int32` are missing, so
  they fail (#480).
- `numpy_dtypes` stores `integer` as well as `double` as float32, so a `uint32` ID above 2^24
  loses precision. Two separate checks add a `storage_warnings` entry when that happens.
- `original_dtype` is set by some paths and not others.
- The mapping is implemented four times: `add_column_to_group`, `add_column_to_group_from_polars`,
  `set_column_with_raw_data` and `csv_to_static.get_column_type`.

On the client:

- `processArrayBuffer` ([DataLoaders.ts](../../src/dataloaders/DataLoaders.ts)) decodes `integer`
  as `Float32Array` and `multitext` as `Uint16Array`.
- A third table, `datastore/datatypes.ts`, said `integer` was `Int32Array` and `multitext`
  `Uint8Array`. `sortWorker` used it, so sorting an `integer` column read float32 bit patterns as
  int32 and put negative values in reverse order. It was removed when that bug was fixed;
  [sortColumns.ts](../../src/datastore/sortColumns.ts) now repeats `DataStore.setColumnData`'s
  rules, which is one more copy for this contract to replace.
- Category values are listed in `datasources.json`, which bloats it for gene-sized vocabularies
  (#643).

## Decision

A column's type has three parts. They are defined once and validated on both sides: a zod schema
in TypeScript and its Python counterpart, ideally generated from one source.

- **Stored type:** a zarr v3 data type (bool, int8 to int64, uint8 to uint64, float32, float64,
  string) at the source's precision. Missing values are part of the stored type, not a reason to
  change type. A fill value marks them only where the type reserves a value that cannot be real
  data, such as NaN for floats. Integer and boolean columns carry a validity mask, as AnnData's
  nullable arrays do, because any integer could be a real value.
- **Encoding:** plain, categorical (integer codes plus a categories array kept in the store beside
  them), or multi-valued.
- **Column kind:** numeric, categorical, unique or multi-valued, which is what charts and filters
  ask about. It is derived from the other two, never stored on its own.

`datatype` stays as a legacy field derived from these, so existing code and h5 projects keep
working, and it is what the upgrade and downgrade scripts translate. Any narrowing, such as
float64 to float32 in the DataStore, happens at load time in one place and is reported. Client
code asks about a column's type only through the existing predicates (`isColumnNumeric`,
`isColumnText`), so the migration changes them and not their callers.

## Consequences

- The four Python encoders become one, used by import, upgrade and downgrade.
- The client picks typed arrays from the contract in one place, instead of in `processArrayBuffer`,
  `DataStore` and `sortColumns.ts` separately.
- Category lists move out of `datasources.json` into the store.
- How the DataStore holds integers that do not fit float32, such as IDs above 2^24, is left to the
  implementation; the contract only requires that narrowing is explicit.
