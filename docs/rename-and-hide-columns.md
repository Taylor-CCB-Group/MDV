# Renaming and hiding columns from a script

Two calls, for use when building MDV projects and views.

```python
project.rename_column(datasource, field, new_name)   # -> True if changed, False if already so
project.soft_delete_column(datasource, field)        # -> True if hidden,  False if already so
```

Both take the column's **`field`** (its identifier), not its display name. `field` never changes,
which is why renaming can't break a chart. Both are safe to re-run.

Hiding sets a tombstone: the column stays in the datasource and in exported files, and MDV stops
showing it. It is reversible.

---

## The one rule: hide columns before you build views

This applies to **hiding only**. `soft_delete_column` refuses while a saved view still references
the column, so hiding has to happen before that view exists.

**`rename_column` is never blocked.** It changes only the display label, and every chart resolves
columns by `field`, which the rename leaves untouched — so you can relabel a column at any point,
including on a project whose views already exist. That is the normal way to fix labels on a
project someone else built.

If you create datasources with `add_datasource`, pass `add_to_view=None` — the default builds a
view listing every column, which blocks every deletion:

```python
project.add_datasource("cells", df, add_to_view=None)   # no view yet
project.rename_column("cells", "pct_counts_mt", "Mito %")
project.soft_delete_column("cells", "scratch_col")
project.create_view_with_all_datasources("default", make_default=True)
```

A project made by `mdvtools convert-scanpy --delete_existing` already starts with an empty view,
so nothing blocks you there.

Views generated afterwards skip hidden columns automatically.

---

## Steps

**1. List the real column names.** Write scripts to a file rather than `python -c`; quoting in
one-liners is a trap.

```bash
cat > /app/mdv/_scripts/list_columns.py <<'PY'
from mdvtools.mdvproject import MDVProject
p = MDVProject("/app/mdv/my_project")
for ds in p.datasources:
    print("==", ds["name"], "-", len(ds["columns"]), "columns")
    for c in ds["columns"]:
        flag = " [hidden]" if c.get("deleted") else ""
        print("  ", c["field"], "|", c["name"], "|", c["datatype"], flag)
PY
```

```bash
cd /app/python && uv run -- python /app/mdv/_scripts/list_columns.py
```

**2. Write the curation script.** A mapping shared across projects will name columns some
projects don't have, so skip those instead of aborting the run.

```python
from mdvtools.mdvproject import MDVProject

RENAME = {
    "n_genes_by_counts": "Gene count",
    "pct_counts_mt":     "Mito %",
}
HIDE = ["scratch_col", "old_clusters"]

p = MDVProject("/app/mdv/my_project")

for field, label in RENAME.items():
    try:
        print("renamed" if p.rename_column("cells", field, label) else "already", field)
    except AttributeError as e:
        print("skipped", field, ":", e)     # not in this project, or not renameable
    except ValueError as e:
        print("FAILED", field, ":", e)      # empty or clashing label - fix the mapping

for field in HIDE:
    try:
        print("hid" if p.soft_delete_column("cells", field) else "already", field)
    except AttributeError as e:
        print("skipped", field, ":", e)
    except ValueError as e:
        print("BLOCKED", field, ":", e)     # a view already uses it - see the rule above
```

The two `except` clauses tell you whether to care. **`AttributeError`** = the column isn't there
or can't be curated — usually fine to skip. **`ValueError`** = you asked for something
contradictory — fix the script.

**3. Re-run step 1** and check `field` is unchanged, `name` is what you wanted, and `[hidden]`
appears only where you meant it.

**4. Build the views.** Reference columns by **`field`**, never by the new display name — chart
configs resolve by field. `create_view_with_all_datasources` handles both automatically.

---

## Over HTTP instead

Get `<id>` from the project's URL in the browser.

```bash
curl -s -X POST http://localhost:5055/project/<id>/rename_column -H 'Content-Type: application/json' -d '{"datasource":"cells","field":"pct_counts_mt","name":"Mito %"}'
```

```bash
curl -s -X POST http://localhost:5055/project/<id>/soft_delete_column -H 'Content-Type: application/json' -d '{"datasource":"cells","field":"scratch_col"}'
```

`200 {"changed":true}` = done · `200 {"changed":false}` = already in that state ·
`400 {"success":false,"error":"..."}` = rejected, with the same message the Python call raises.

---

## Common errors

| Message | Meaning |
|---|---|
| `column X not found in Y datasource` | wrong field — you may be passing the display name |
| `column X already exists in Y datasource` | another visible column already uses that label |
| `Column name is required` | empty label |
| `column X ... is a spatial (sgindex) column` | geometry column, not curatable |
| `deletion of column X ... is blocked: it is still used by other charts or views (view: chart)` | hide the column before that view exists, or edit the chart first |

---

## Worth knowing

- **Hiding doesn't remove data.** `get_datasource_as_dataframe()` and exported files still
  include hidden columns. A script that counts columns will count them.
- **The dataframe stays keyed by `field`.** After renaming, it's still `df["n_genes_by_counts"]`.
- **The block can over-report.** The check looks for the field string anywhere in a chart config,
  so a chart *titled* the same as a field will block that deletion. Deliberate: a false "blocked"
  is recoverable, a false "clear" leaves a chart pointing at a missing column.
- **Run curation when nobody has the project open.** A browser tab loaded beforehand can save
  stale metadata back over your change.
- **Don't hide structural columns** — a datasource's `links` `name_column`, an image set's
  `key_column`, or embedding columns a scatter plot needs. The block only inspects chart configs,
  so it won't stop you.

## Undo

`set_column_metadata` bypasses every check — right for repair, wrong for curation.

```python
project.set_column_metadata("cells", "n_genes_by_counts", "name", "n_genes_by_counts")  # undo rename
project.set_column_metadata("cells", "scratch_col", "deleted", False)                    # un-hide
```
