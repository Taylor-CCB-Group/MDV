# Spatial Connectivity Map - Metadata-Driven Network Chart

## Overview

The **Spatial Connectivity Map** (CellNetworkChart) creates interactive force-directed network visualizations with Python-configured metadata. This chart requires setup but provides advanced features like ego networks and spatial integration.

**Chart Type**: `cell_network_chart` (displays as "Spatial Connectivity Map" in GUI)  
**Status**: ✅ Available in MDV GUI **only after Python setup**  
**Location**: `/app/src/charts/CellNetworkChart.js`

> **⚠️ Important**: This chart **will not appear** in the Add Chart dialog unless you configure the `interactions` metadata using Python (see setup below). For quick exploration without Python setup, use the [Network Graph](/app/docs/NETWORK_GRAPH_GUIDE.md) instead!

---

## Why Not Visible in GUI?

**The "Spatial Connectivity Map" chart requires Python configuration first.** It will **not appear** in the Add Chart dialog until you run `setup_ligand_network()` in Python.

This is by design - the chart requires specific metadata that must be configured programmatically. If you want a chart that works immediately without setup, use [Network Graph](NETWORK_GRAPH_GUIDE.md) instead.

---

## Quick Start (5 minutes)

### Step 1: Prepare Your Data

Your interaction data needs these columns:

| Required | Column | Type | Example Values |
|----------|--------|------|----------------|
| ✅ | Pivot/Category | Text | "VEGFA", "TGFB1", "immune_cancer" |
| ✅ | Source Cell ID | Text | "T_cell_1", "Cancer_A" |
| ✅ | Target Cell ID | Text | "Endothelial_1", "Fibroblast_2" |
| ✅ | Interaction Score | Number | 0.85, 0.72, 0.91 |
| ⭕ | Spatial Distance | Number | 12.5, 8.3, 15.2 (µm) |
| ⭕ | P-value/FDR | Number | 0.001, 0.003, 0.0001 |
| ⭕ | Cell Count | Number | 15, 22, 8 |
| ⭕ | Cell Type | Text | "T_cell", "Cancer", "Fibroblast" |

### Step 2: Python Setup (2 lines of code!)

```python
from mdvtools import MDVProject
from mdvtools.network_helpers import setup_ligand_network
import pandas as pd

# Your interaction data
df = pd.read_csv("interactions.csv")  # or from CellPhoneDB, etc.

# Create project
project = MDVProject("my_network_analysis")
project.add_datasource("interactions", data=df, size=len(df))

# Configure network (automatically detects best columns)
setup_ligand_network(
    project,
    interaction_datasource="interactions",
    ligand_column="ligand_type",           # or communication_type, pathway, etc.
    source_cell_column="source_cell_id",    # sender/from cell
    target_cell_column="target_cell_id",    # receiver/to cell
    interaction_score_column="score",       # strength/probability
    spatial_distance_column="distance",     # optional
    pvalue_column="pvalue"                  # optional
)

project.save()
```

### Step 3: Open in GUI

1. Launch MDV and load your project
2. View should show "Ligand Network Analysis" (or your custom name)
3. **Click "Add Chart" button** (top of interface)
4. **Select datasource**: Your interaction datasource (e.g., "interactions")
5. **Chart type**: Select **"Spatial Connectivity Map"** from dropdown
6. **Filter dropdown**: Select which category to visualize (e.g., "VEGFA")
7. **Click "Add Chart"** → Network appears! 🎉

---

## Visual Encoding

The network uses 4 visual channels to encode your data:

| Visual Property | Maps To | Default Column | Adjustable? |
|----------------|---------|----------------|-------------|
| **Link Thickness** | Interaction strength | `interaction_score` | ✅ Yes |
| **Link Length** | Spatial distance (or score) | `spatial_distance` | ✅ Yes |
| **Link Color** | Statistical significance | `pvalue` (or score) | ✅ Yes |
| **Node Size** | Cell count / interactions | `cell_count` (or score) | ✅ Yes |
| **Node Color** | Cell type / state | `cell_type` (optional) | ✅ Yes |

---

## Interactive Features

### Mouse Interactions
- **Drag nodes**: Manually reposition cells
- **Click link**: Highlights interaction (propagates to other charts)
- **Click node**: Highlights all interactions for that cell
- **Scroll wheel**: Zoom in/out (if enabled)

### Settings Panel (⚙️ icon on chart)

**Physics Controls:**
- **Link Strength** (0-2): How tightly links pull nodes together
  - Low (0.3): Loose, spread out
  - High (1.5): Tight clusters
- **Node Repulsion** (-1000 to 0): How strongly nodes push apart
  - Strong (-800): Nodes far apart
  - Weak (-100): Nodes can overlap

**Visual Controls:**
- **Link Thickness**: Change column, domain (min/max), range (pixels)
- **Link Length**: Change column, domain, range
- **Link Color**: Choose color scheme (red, blue-yellow, etc.), adjust domain
- **Node Size**: Adjust radius range (1-50px)
- **Node Color**: By cells / by state / by type
- **Show Directionality**: Toggle arrows on links

**Layout Controls:**
- **Center Cell**: Select a cell to focus on (ego network mode)
- **Number of Levels**: Show 1 or 2 hops from center
- **Label Size**: Adjust node label font size (5-20px)

---

## Advanced Usage

### 1. Ego Network Mode

Focus on one cell and its neighborhood:

```javascript
// In chart settings:
Center Cell: "T_cell_1"
Number of Levels: 2
```

Shows:
- **Level 1**: `T_cell_1` and all directly connected cells
- **Level 2**: Cells connected to Level 1 cells

Great for exploring local neighborhoods in large networks!

### 2. Multiple Network Views

Create separate views for each ligand/pathway:

```python
ligand_types = df['ligand_type'].unique()

for ligand in ligand_types:
    # Each ligand gets its own view
    view_name = f"{ligand} Signaling Network"
    # ... configure view
```

### 3. Link to Spatial Data

If you have cell coordinates, link network to spatial view:

```python
# Load spatial coordinates
cells_df = pd.read_csv("cell_positions.csv")  # must have x, y columns
project.add_datasource("cells", data=cells_df, size=len(cells_df))

# Set up spatial regions
project.set_region_data(
    "cells",
    {"roi_1": {"width": 1000, "height": 1000}},
    region_field="roi",
    default_color="cell_type",
    position_fields=["x", "y"],
    scale_unit="µm",
    scale=1.0
)

# Link to network
setup_ligand_network(
    project,
    interaction_datasource="interactions",
    cells_datasource="cells",  # ← connects network to spatial
    ...
)
```

**Result**: Click network link → Highlights cells in spatial viewer!

### 4. Bidirectional Link Handling

Network automatically handles bidirectional interactions:

```
If both exist:
  A → B (score: 0.85)
  B → A (score: 0.45)

Chart keeps only: A → B (stronger)
```

This prevents:
- Visual clutter from overlapping lines
- Double-counting of interactions
- Confusion about directionality

---

## Integration with Analysis Tools

### CellPhoneDB Output

```python
import pandas as pd

# Load CellPhoneDB results
cpdb_means = pd.read_csv("means.txt", sep="\t")
cpdb_pvals = pd.read_csv("pvalues.txt", sep="\t")

# Transform to MDV format
interactions = []
for idx, row in cpdb_means.iterrows():
    ligand_receptor = row['interacting_pair']
    ligand, receptor = ligand_receptor.split('_')
    
    # Each cell type pair is a separate row
    for col in cpdb_means.columns[11:]:  # Skip metadata columns
        cell_a, cell_b = col.split('|')
        interactions.append({
            'ligand_type': ligand,
            'receptor_type': receptor,
            'source_cell_id': cell_a,
            'target_cell_id': cell_b,
            'interaction_score': row[col],
            'pvalue': cpdb_pvals.loc[idx, col]
        })

df = pd.DataFrame(interactions)
# Filter significant interactions
df = df[df['pvalue'] < 0.05]

# Load into MDV
project.add_datasource("cellphonedb", data=df, size=len(df))
setup_ligand_network(project, "cellphonedb", ligand_column="ligand_type", ...)
```

### NicheNet / LIANA / CellChat

Similar transformation - map their output columns to MDV format:

```python
# Generic mapping function
def map_communication_results(df, tool="nicheNet"):
    """Map different tool outputs to MDV format"""
    
    column_mappings = {
        "nicheNet": {
            'ligand': 'ligand_type',
            'receptor': 'receptor_type',
            'sender': 'source_cell_id',
            'receiver': 'target_cell_id',
            'activity': 'interaction_score'
        },
        "LIANA": {
            'ligand_complex': 'ligand_type',
            'receptor_complex': 'receptor_type',
            'source': 'source_cell_id',
            'target': 'target_cell_id',
            'aggregate_rank': 'interaction_score'
        },
        # Add more tools...
    }
    
    mapping = column_mappings[tool]
    return df.rename(columns=mapping)
```

---

## Troubleshooting

### Chart doesn't appear in Add Chart dialog

**Problem**: "Spatial Connectivity Map" not in chart type dropdown

**Solutions**:
1. ✅ Check datasource has `interactions` metadata:
   ```python
   md = project.get_datasource_metadata("your_datasource")
   print(md.get("interactions"))  # Should NOT be None
   ```
2. ✅ Run `setup_ligand_network()` to configure metadata
3. ✅ Reload project in GUI (close & reopen)

### Network shows no nodes/links

**Problem**: Chart loads but is empty

**Solutions**:
1. ✅ Check you selected a pivot category (ligand type) in Add Chart dialog
2. ✅ Verify data has rows for that category:
   ```python
   df[df['ligand_type'] == 'VEGFA']  # Should have rows
   ```
3. ✅ Check filters aren't hiding everything (clear filters in GUI)
4. ✅ Verify source_cell != target_cell (no self-loops)

### Links are too long/short

**Problem**: Network is too spread out or too compact

**Solutions**:
1. ✅ Adjust **Link Strength** in settings (lower = longer links)
2. ✅ Adjust **Node Repulsion** (more negative = more spread)
3. ✅ Change **Link Length Domain** in settings
4. ✅ Use different column for link length (e.g., constant value)

### Can't see node labels

**Problem**: Labels are missing or too small

**Solutions**:
1. ✅ Increase **Label Size** in settings (5-20px)
2. ✅ Zoom in on chart
3. ✅ Make chart larger (resize chart panel)

### Links all same color

**Problem**: Link color doesn't vary

**Solutions**:
1. ✅ Check p-value/color column has variation:
   ```python
   df['pvalue'].describe()  # Should have range
   ```
2. ✅ In settings: **Link Color Domain** → adjust min/max
3. ✅ Try different color scheme (red, blue-yellow, etc.)

---

## Example Data Formats

### Minimal Example (3 required columns)
```csv
interaction_type,source_cell,target_cell,score
immune_signal,T_cell_1,Cancer_A,0.85
immune_signal,T_cell_2,Cancer_A,0.72
stromal_signal,Fibroblast_1,Cancer_B,0.91
```

### Full Example (all optional columns)
```csv
ligand_type,receptor_type,source_cell_id,target_cell_id,source_type,target_type,score,distance,pvalue,fdr,count
VEGFA,KDR,T_cell_1,Endothelial_1,T_cell,Endothelial,0.85,12.5,0.001,0.005,15
VEGFA,FLT1,T_cell_1,Endothelial_2,T_cell,Endothelial,0.72,8.3,0.003,0.012,22
TGFB1,TGFBR1,Cancer_A,Fibroblast_1,Cancer,Fibroblast,0.91,15.2,0.0001,0.001,8
```

---

## Complete Python API

```python
from mdvtools.network_helpers import (
    setup_ligand_network,           # Main configuration function
    create_ligand_network_view,     # Create view with network
    setup_cell_communication_network # Alternative naming
)

# Full parameter list
setup_ligand_network(
    project,                              # MDVProject instance
    interaction_datasource="interactions", # Name of datasource
    cells_datasource="cells",             # Optional: for spatial linking
    
    # Column mappings (required)
    ligand_column="ligand_type",          # Category to filter by
    source_cell_column="source_cell_id",  # Source nodes
    target_cell_column="target_cell_id",  # Target nodes
    interaction_score_column="score",     # Numeric strength
    
    # Visual encoding (optional)
    spatial_distance_column="distance",   # For link length
    pvalue_column="pvalue",               # For link color
    cell_count_column="count",            # For node size
    cell_type_column="cell_type",         # For node color
    
    # View creation
    create_view=True,                     # Auto-create view
    view_name="Ligand Network Analysis"   # View name
)
```

---

## Tips & Best Practices

### 1. **Start with one ligand/pathway**
- Don't try to visualize all interactions at once
- Use pivot column to filter to specific types
- Create separate views for different pathways

### 2. **Adjust physics first**
- Get the layout right before tweaking colors/sizes
- Link strength: 0.3-0.5 for large networks, 1.0+ for small
- Node repulsion: -500 for medium networks

### 3. **Use meaningful column names**
- GUI shows column names in legends/tooltips
- Rename columns before loading for clarity
- Example: "pvalue" → "Statistical Significance"

### 4. **Filter before visualizing**
- Pre-filter to significant interactions (p < 0.05)
- Remove weak interactions (score < threshold)
- Reduces clutter, improves performance

### 5. **Link to spatial when possible**
- Provides crucial context about cell locations
- Click network → see cells highlighted in space
- Reveals spatial patterns in signaling

---

## Need More Help?

- **Full Example**: See `/app/python/mdvtools/examples/ligand_network_example.py`
- **Source Code**: `/app/src/charts/CellNetworkChart.js`
- **Python API**: `/app/python/mdvtools/network_helpers.py`
- **MDV Docs**: See main MDV documentation

---

## Summary

**The CellNetworkChart (Spatial Connectivity Map) is already available in MDV!**

Just:
1. ✅ Prepare interaction data (4 required columns)
2. ✅ Run `setup_ligand_network()` in Python
3. ✅ Open GUI → Add Chart → Select "Spatial Connectivity Map"
4. ✅ Choose ligand/pathway type
5. ✅ Explore your network!

Happy networking! 🕸️

