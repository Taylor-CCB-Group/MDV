# Flexible Network Chart - Zero-Friction Network Visualization

## Overview

The **Network Graph (Flexible)** chart works with **ANY table** - no Python setup, no metadata configuration required. Just select your columns in the GUI and visualize!

**Chart Type**: `flexible_network_chart`  
**Display Name**: "Network Graph (Flexible)"  
**Status**: ✅ Available now

---

## Quick Start (2 minutes, GUI only!)

### Your Data (any CSV/table):
```csv
from_cell,to_cell,interaction_strength,ligand,distance,pvalue,cell_count,cell_type
T_cell_1,Cancer_A,0.85,VEGFA,12.5,0.001,15,T_cell
T_cell_2,Cancer_B,0.72,VEGFA,8.3,0.003,22,T_cell
Cancer_A,Fibroblast_1,0.91,TGFB1,15.2,0.0001,8,Cancer
```

### Steps:
1. **Load data** in MDV (via Python or file upload)
2. **Click "Add Chart"** button
3. **Select your datasource** (any table with network data)
4. **Chart type**: Select **"Network Graph (Flexible)"**
5. **Select columns** via dropdowns:
   - Source Node → `from_cell`
   - Target Node → `to_cell`
   - Edge Weight → `interaction_strength`
   - Category Filter (optional) → `ligand`
   - Edge Length (optional) → `distance`
   - Edge Color (optional) → `pvalue`
   - Node Size (optional) → `cell_count`
   - Node Type (optional) → `cell_type`
6. **Click "Add Chart"** → Network appears!

**That's it! No Python configuration needed!** 🎉

---

## Column Selection Guide

| Column | Required? | Type | Purpose | Example |
|--------|-----------|------|---------|---------|
| **Source Node** | ✅ Required | Text | Start of edge | `cell_1`, `gene_A` |
| **Target Node** | ✅ Required | Text | End of edge | `cell_2`, `gene_B` |
| **Edge Weight** | ✅ Required | Number | Link thickness | `0.85`, `42` |
| **Category Filter** | ⭕ Optional | Text | Filter by type | `VEGFA`, `pathway_1` |
| **Edge Length** | ⭕ Optional | Number | Physical spacing | `12.5` (µm) |
| **Edge Color** | ⭕ Optional | Number | Significance | `0.001` (p-value) |
| **Node Size** | ⭕ Optional | Number | Node radius | `15` (count) |
| **Node Type** | ⭕ Optional | Text | Node coloring | `T_cell`, `Cancer` |

### Minimum Requirements
You only need **3 columns**:
- Source Node (text)
- Target Node (text) 
- Edge Weight (number)

Everything else is optional!

---

## Use Cases

### 1. Ligand-Receptor Interactions (No Python!)
```csv
ligand_cell,receptor_cell,score,ligand_type,distance,pvalue
T_cell_1,Endothelial_1,0.85,VEGFA,12.5,0.001
T_cell_2,Endothelial_2,0.72,VEGFA,8.3,0.003
```

**In GUI:**
- Source: `ligand_cell`
- Target: `receptor_cell`
- Weight: `score`
- Category: `ligand_type`
- Edge Color: `pvalue`

Done! Filter by VEGFA in settings to see only VEGFA interactions.

### 2. Gene Regulatory Networks
```csv
transcription_factor,target_gene,regulatory_score,regulation_type
STAT3,IL6,0.92,activation
STAT3,SOCS3,0.87,activation
TP53,MDM2,0.95,activation
```

**In GUI:**
- Source: `transcription_factor`
- Target: `target_gene`
- Weight: `regulatory_score`
- Category: `regulation_type`

### 3. Cell-Cell Communication
```csv
sender,receiver,communication_score,pathway,cell_distance
Neuron_1,Astrocyte_1,0.88,glutamate_signaling,5.2
Neuron_2,Microglia_1,0.76,cytokine_signaling,8.7
```

**In GUI:**
- Source: `sender`
- Target: `receiver`
- Weight: `communication_score`
- Category: `pathway`
- Edge Length: `cell_distance`

### 4. Social Networks / Any Graph
```csv
user_from,user_to,connection_strength,relationship_type
alice,bob,0.95,friend
bob,charlie,0.82,colleague
charlie,alice,0.90,friend
```

Works with **literally any network data**!

---

## Visual Encoding

| Visual Property | Maps To | Selected Column |
|----------------|---------|-----------------|
| **Link Thickness** | Edge strength | Edge Weight (always) |
| **Link Length** | Physical distance | Edge Length (if provided) |
| **Link Color** | Significance/value | Edge Color (if provided) |
| **Node Size** | Count/importance | Node Size (if provided) |
| **Node Color** | Type/category | Node Type (if provided) |

---

## Interactive Features

### Mouse Controls
- **Drag nodes**: Reposition manually
- **Click link**: Highlight edge (shows in other charts)
- **Click node**: Highlight all edges for that node
- **Hover link**: Shows tooltip with source → target and weight

### Settings Panel (⚙️)

**Physics:**
- **Link Strength** (0-2): How tightly connected
  - 0.3 = Loose, spread out
  - 1.5 = Tight, compact
- **Node Repulsion** (-1000 to 0): How much nodes push apart
  - -800 = Very spread out
  - -100 = Nodes can get close

**Display:**
- **Node Radius** (3-30): Size of nodes (if Node Size not selected)
- **Show Node Labels**: Toggle text labels
- **Show Directionality**: Toggle arrows on edges

---

## Comparison: Flexible vs. Metadata-Driven

| Feature | Flexible Network | Spatial Connectivity Map |
|---------|-----------------|-------------------------|
| **Setup** | None - just GUI | Python metadata required |
| **Works with** | ANY table | Only configured datasources |
| **Column selection** | GUI dropdowns (8 options) | Pre-configured in Python |
| **Friction** | Zero | Medium |
| **Flexibility** | Maximum | Limited |
| **Use case** | Exploration, varied data | Production, standardized data |

**Use Flexible Network when:**
- ✅ You want to quickly explore data
- ✅ Your data structure varies
- ✅ You don't want to write Python code
- ✅ You're working with non-standard formats

**Use Spatial Connectivity Map when:**
- ✅ You have standardized data pipelines
- ✅ You want validated, consistent configurations
- ✅ You're building production workflows

---

## Advanced: Category Filtering

If you select a **Category Filter** column (e.g., `ligand_type`):

1. Add the chart normally
2. Open **Settings Panel** (⚙️)
3. Find chart config (may need to edit JSON for now)
4. Set `category_filter` to filter value:

```javascript
{
  "type": "flexible_network_chart",
  "param": [...],
  "category_filter": "VEGFA"  // Only show VEGFA interactions
}
```

Or filter the data before loading into MDV!

---

## Tips & Best Practices

### 1. Start Simple
- Begin with just 3 required columns
- Add optional columns one at a time
- See what works best for your data

### 2. Name Columns Clearly
- Column names appear in tooltips and legends
- Use descriptive names: "P-value" not "col5"
- Rename before loading for clarity

### 3. Pre-filter Large Networks
- Networks with >1000 edges can be slow
- Filter to significant interactions first
- Use Category Filter to show subsets

### 4. Adjust Physics First
- Get layout right before colors/sizes
- Link Strength 0.3-0.7 for most networks
- Node Repulsion -400 is a good starting point

### 5. Use Edge Color for P-values
- Lower p-value = darker color
- Quickly identify significant interactions
- Set range in settings if needed

---

## Example Workflows

### From CellPhoneDB Output

```python
# Load CellPhoneDB results
import pandas as pd
from mdvtools import MDVProject

means = pd.read_csv("means.txt", sep="\t")

# Reshape for network format
interactions = []
for idx, row in means.iterrows():
    ligand_receptor = row['interacting_pair']
    for col in means.columns[11:]:
        cell_pair = col.split('|')
        interactions.append({
            'source_cell': cell_pair[0],
            'target_cell': cell_pair[1],
            'interaction_strength': row[col],
            'ligand_receptor_pair': ligand_receptor
        })

df = pd.DataFrame(interactions)
df = df[df['interaction_strength'] > 0]  # Remove zeros

# Load into MDV
project = MDVProject("cellphonedb_network")
project.add_datasource("interactions", data=df, size=len(df))
project.save()
```

Then in GUI:
- Source: `source_cell`
- Target: `target_cell`  
- Weight: `interaction_strength`
- Category: `ligand_receptor_pair`

### From Simple CSV

```python
# Literally just load your CSV
import pandas as pd
from mdvtools import MDVProject

df = pd.read_csv("my_network.csv")

project = MDVProject("my_network_viz")
project.add_datasource("network", data=df, size=len(df))
project.save()
```

Open in GUI → Add Chart → Select columns → Done!

---

## Troubleshooting

### Network is empty
- ✅ Check source ≠ target (no self-loops removed automatically)
- ✅ Verify data has rows after filtering
- ✅ Check column selections are correct

### Network is too spread out
- ✅ Increase Link Strength (0.5 → 1.5)
- ✅ Decrease Node Repulsion (-400 → -200)
- ✅ Make chart panel smaller (forces compact layout)

### Network is too compact
- ✅ Decrease Link Strength (0.5 → 0.2)
- ✅ Increase Node Repulsion (-400 → -800)
- ✅ Use Edge Length column with larger values

### Can't see edges
- ✅ Increase Edge Weight range in settings
- ✅ Check weight values aren't all very small
- ✅ Zoom out or make chart larger

### Nodes all same size
- ✅ Select Node Size column in params
- ✅ If already selected, check column has variation
- ✅ Adjust node size range in settings

---

## Complete Example

### Data: `interactions.csv`
```csv
from,to,score,type,distance,pval,count,category
A,B,0.8,ligand,10,0.01,5,immune
A,C,0.9,receptor,15,0.001,8,immune
B,C,0.7,ligand,12,0.05,3,immune
D,E,0.95,receptor,8,0.0001,12,stromal
```

### Python (Optional - just to load):
```python
from mdvtools import MDVProject
import pandas as pd

df = pd.read_csv("interactions.csv")
project = MDVProject("flex_network")
project.add_datasource("data", data=df, size=len(df))
project.save()
```

### GUI:
1. Open MDV, load "flex_network"
2. Add Chart → datasource: "data"
3. Type: "Network Graph (Flexible)"
4. Columns:
   - Source: `from`
   - Target: `to`
   - Weight: `score`
   - Category: `type`
   - Edge Length: `distance`
   - Edge Color: `pval`
   - Node Size: `count`
   - Node Type: `category`
5. Add Chart!

**Result**: Interactive network with:
- Thickness = score
- Length = distance
- Color = pval (red=significant)
- Node size = count
- Node color = category (immune vs stromal)

---

## Summary

**The Flexible Network Chart removes ALL friction!**

✅ **No Python setup** - just GUI selections  
✅ **Works with ANY table** - no special format  
✅ **8 column dropdowns** - pick what you need  
✅ **3 required, 5 optional** - start simple  
✅ **Instant visualization** - no metadata  

**Perfect for:**
- Quick exploration
- Varied data formats
- Non-standard analyses
- Users without Python expertise
- Rapid prototyping

**Go from CSV to network in 2 minutes!** 🚀

