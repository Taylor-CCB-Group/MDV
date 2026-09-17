# MDVTools Examples

This directory contains example scripts demonstrating various MDV features.

## Network Visualization

### `ligand_network_example.py`

Complete example showing how to create interactive ligand-receptor network visualizations.

**Run it:**
```bash
python ligand_network_example.py
```

**Creates:**
- Example ligand-receptor interaction data
- MDV project configured for network visualization  
- View with network chart ready to use

**Then:**
1. Open MDV GUI
2. Load the created project
3. Add "Spatial Connectivity Map" chart
4. Select a ligand type (VEGFA, TGFB1, etc.)
5. Explore the interactive network!

**Features demonstrated:**
- Network configuration with `setup_ligand_network()`
- Visual encoding (link thickness, color, node size)
- Linking to spatial data
- View creation with instructional text

## Column Curation

### `rename_delete_columns_example.py`

Renaming and hiding columns from a script, and the ordering rule that makes it work.

**Run it:**
```bash
python rename_delete_columns_example.py
```

**Creates:**
- A small `cells` datasource with the sort of column names a conversion produces
- Four renamed display labels and two hidden columns
- A default view built *after* curation, listing only the surviving columns

**Features demonstrated:**
- `rename_column()` and `soft_delete_column()` - return values and both exception types
- `add_datasource(..., add_to_view=None)` - why curation must happen before any view exists
- `create_view_with_all_datasources()` skipping hidden columns automatically
- The view-reference guard refusing a deletion once a chart uses the column

See `docs/rename-and-hide-columns.md` for the full guide.

## More Examples Coming Soon

- Spatial data integration
- Multi-omics visualization
- Custom chart configurations
- Advanced filtering and selection

