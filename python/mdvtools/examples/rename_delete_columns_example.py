"""
Renaming and hiding columns from a script.

Two operations on column metadata:

  - rename_column(datasource, field, new_name)  - change the display label
  - soft_delete_column(datasource, field)       - hide the column from MDV

Both address the column by its `field` (its identifier), never by its display name.
`field` never changes, which is why renaming cannot break a chart. Hiding sets a
tombstone: the column stays in the datasource and in exported files, and MDV stops
showing it. Both are reversible and safe to re-run.

THE ORDERING RULE, which this example exists to demonstrate - it applies to
HIDING only:

    soft_delete_column refuses while a saved view still references the column.
    So hide columns BEFORE building views.

    add_datasource defaults to add_to_view="default", which generates a table plot
    listing every column - and that view would then block every deletion. Passing
    add_to_view=None is what makes hiding possible.

    rename_column is never blocked. It changes only the display label, and charts
    resolve columns by `field`, which a rename leaves untouched - so labels can be
    fixed at any time, including on a project whose views already exist.

See docs/rename-and-hide-columns.md for the full guide.
"""

import pandas as pd
import numpy as np
from mdvtools import MDVProject

PROJECT = "rename_delete_columns_example"

# =============================================================================
# 1. Some data with the sort of column names a conversion tends to produce
# =============================================================================

np.random.seed(0)
n = 200

cells = pd.DataFrame({
    "n_genes_by_counts": np.random.randint(500, 4000, n),
    "total_counts":      np.random.randint(1000, 20000, n),
    "pct_counts_mt":     np.round(np.random.uniform(0, 0.15, n), 3),
    "leiden":            np.random.choice(["0", "1", "2", "3"], n),
    # working columns nobody outside the analysis needs to see
    "scratch_qc":        np.random.uniform(0, 1, n),
    "old_clusters":      np.random.choice(["a", "b"], n),
})

# =============================================================================
# 2. Create the project WITHOUT a view - this is the step that matters
# =============================================================================

project = MDVProject(PROJECT, delete_existing=True)
project.add_datasource("cells", cells, add_to_view=None)


def show(label):
    """Print each column as: field -> display name [hidden]"""
    columns = project.get_datasource_metadata("cells")["columns"]
    print("\n" + label)
    for c in columns:
        flag = "  [hidden]" if c.get("deleted") else ""
        print("   ", c["field"], "->", c["name"], flag)


show("before curation:")

# =============================================================================
# 3. Rename - give the columns labels a reader can interpret
# =============================================================================

RENAME = {
    "n_genes_by_counts": "Gene count",
    "total_counts":      "Total counts",
    "pct_counts_mt":     "Mito %",
    "leiden":            "Cluster",
}

print("\nrenaming:")
for field, label in RENAME.items():
    try:
        changed = project.rename_column("cells", field, label)
        print("   ", "renamed" if changed else "already", field, "->", label)
    except AttributeError as e:
        # not in this project, or not curatable (spatial column, no `field`, ...)
        # a mapping shared across projects will hit this - skip rather than abort
        print("    skipped", field, ":", e)
    except ValueError as e:
        # empty label, or another visible column already uses it - fix the mapping
        print("    FAILED", field, ":", e)

# =============================================================================
# 4. Hide - keep working columns out of the user's way
# =============================================================================

HIDE = ["scratch_qc", "old_clusters"]

print("\nhiding:")
for field in HIDE:
    try:
        changed = project.soft_delete_column("cells", field)
        print("   ", "hid" if changed else "already hidden", field)
    except AttributeError as e:
        print("    skipped", field, ":", e)
    except ValueError as e:
        # a saved view already references it - see the ordering rule above
        print("    BLOCKED", field, ":", e)

show("after curation:")

# =============================================================================
# 5. NOW build the view. Generated views skip hidden columns automatically.
# =============================================================================

view_name = project.create_view_with_all_datasources("default", make_default=True)
project.set_editable(True)

charts = project.get_view(view_name)["initialCharts"]["cells"]
print("\nview '" + view_name + "' created with", len(charts), "chart(s)")

# =============================================================================
# 6. The guard, demonstrated: now that a chart uses these columns, they are
#    protected. This is why curation has to come first.
# =============================================================================

print("\ntrying to hide a column the new view uses:")
try:
    project.soft_delete_column("cells", "pct_counts_mt")
    print("    unexpectedly succeeded - the view should have blocked this")
except ValueError as e:
    print("    blocked, as expected:")
    print("   ", e)

print("\n" + "=" * 60)
print("Project written to:", project.dir)
print("=" * 60)
print("\nNotes:")
print("  - `field` never changed - charts, links and column groups still resolve")
print("  - hidden columns are still in the data and in exported files")
print("  - to un-hide: set_column_metadata('cells', <field>, 'deleted', False)")
