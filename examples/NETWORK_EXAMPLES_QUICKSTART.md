# Network Chart Examples - Quick Start 🚀

## 🎯 Instant Demo

**Want to see the network chart in action right now?**

```bash
cd /app/examples
python load_network_examples.py
```

Then **open MDV GUI** → Load project "network_examples" → Follow instructions in each view!

---

## 📁 Example CSV Files

Located in `/app/examples/data/`:

### 1. **ligand_receptor_network_example.csv** ⭐ RECOMMENDED

**Full-featured demo** showing ALL chart capabilities.

```csv
source_cell,target_cell,interaction_score,ligand_type,spatial_distance,pvalue,cell_count,cell_type
T_cell_1,Endothelial_1,0.89,VEGFA,12.5,0.0001,18,T_cell
T_cell_1,Endothelial_2,0.76,VEGFA,8.3,0.0023,18,T_cell
...
```

**What it shows:**
- 30 interactions between 20 cells
- 7 ligand types (VEGFA, TGFB1, PDGFB, etc.)
- 5 cell types (T_cell, Cancer, Fibroblast, etc.)
- Spatial distances, p-values, cell counts

**All 8 columns mapped to visual properties!**

---

### 2. **simple_network_example.csv**

**Minimal demo** with only 3 required columns.

```csv
from_node,to_node,weight
A,B,0.85
A,C,0.92
B,C,0.78
...
```

**What it shows:**
- 11 edges, 7 nodes (A through G)
- Just the basics: source, target, weight
- Perfect for understanding the minimum requirements

---

### 3. **gene_regulatory_network_example.csv**

**Gene regulation demo** for biological networks.

```csv
transcription_factor,target_gene,regulatory_score,regulation_type,binding_affinity
TP53,MDM2,0.95,activation,8.2
TP53,BAX,0.88,activation,7.5
...
```

**What it shows:**
- 18 TF-gene interactions
- 5 transcription factors (TP53, STAT3, NFKB, MYC, HIF1A)
- Regulatory relationships with binding scores

---

## 🎬 Quick Demo (2 minutes)

### Option A: Automatic (Easiest!)

```bash
# Run this script
cd /app/examples
python load_network_examples.py

# Then open MDV GUI and load "network_examples" project
# Follow instructions in each view!
```

### Option B: Manual

1. **Load CSV in Python:**
```python
from mdvtools import MDVProject
import pandas as pd

# Pick any CSV
df = pd.read_csv("/app/examples/data/simple_network_example.csv")

project = MDVProject("my_demo")
project.add_datasource("network", data=df, size=len(df))
project.save()
```

2. **Open in GUI:**
   - Launch MDV
   - Load "my_demo" project
   - Click "Add Chart"
   - Select "network" datasource
   - Chart type: **"Network Graph (Flexible)"**
   - Select columns from dropdowns
   - Click "Add Chart"

3. **See your network!** 🎉

---

## 📊 Column Selection Guide

### For `ligand_receptor_network_example.csv` (Full Example)

In the Add Chart dialog, select:

| Dropdown | Column to Select |
|----------|------------------|
| Source Node | `source_cell` |
| Target Node | `target_cell` |
| Edge Weight | `interaction_score` |
| Category Filter | `ligand_type` |
| Edge Length | `spatial_distance` |
| Edge Color | `pvalue` |
| Node Size | `cell_count` |
| Node Type | `cell_type` |

**Result:** Fully-featured network with all visual encodings!

---

### For `simple_network_example.csv` (Minimal)

In the Add Chart dialog, select:

| Dropdown | Column to Select |
|----------|------------------|
| Source Node | `from_node` |
| Target Node | `to_node` |
| Edge Weight | `weight` |
| *(Leave others empty)* | |

**Result:** Clean, simple network!

---

### For `gene_regulatory_network_example.csv`

In the Add Chart dialog, select:

| Dropdown | Column to Select |
|----------|------------------|
| Source Node | `transcription_factor` |
| Target Node | `target_gene` |
| Edge Weight | `regulatory_score` |
| Category Filter | `regulation_type` *(optional)* |
| Edge Color | `binding_affinity` *(optional)* |

**Result:** Gene regulatory network!

**Tip:** After creating, enable "Show Directionality" in settings to see TF → gene arrows.

---

## 🎨 What the Visuals Mean

| Visual Property | Maps To | Example |
|----------------|---------|---------|
| **Link thickness** | Edge weight/score | Thick = strong interaction |
| **Link length** | Spatial distance | Short = cells close together |
| **Link color** | P-value/significance | Red = significant, gray = not |
| **Node size** | Cell count/importance | Large = many interactions |
| **Node color** | Cell/node type | Different colors = different types |

---

## 🔧 Interactive Features

### After Creating a Chart:

**Drag nodes** - Click and drag to rearrange manually

**Click edges** - Highlights that interaction (shows in other charts too!)

**Click nodes** - Highlights all edges for that node

**Settings (⚙️ icon):**
- Link Strength (0-2): How tightly connected - try 0.5
- Node Repulsion (-1000 to 0): How spread out - try -400
- Show Labels: Toggle node text
- Show Directionality: Toggle arrows

**Filters:**
- Use MDV's data filters to show subsets
- Filter by `ligand_type` to see specific pathways
- Filter by `cell_type` to focus on interactions

---

## 💡 Tips for Best Results

### 1. Start with Simple Example
- Understand basics with just 3 columns
- Add complexity gradually

### 2. Adjust Physics
- **Spread out**: Decrease Link Strength (0.3), increase Node Repulsion (-600)
- **Compact**: Increase Link Strength (1.0), decrease Node Repulsion (-200)

### 3. Use Colors Wisely
- Edge color for p-values: Red = significant
- Node color for cell types: Different categories
- Helps identify patterns quickly

### 4. Filter Your Data
- Large networks (>50 edges) can be cluttered
- Filter to significant interactions first
- Use Category Filter column to show subsets

---

## 📝 Creating Your Own CSV

### Minimum (3 columns):
```csv
from,to,weight
nodeA,nodeB,0.85
nodeB,nodeC,0.72
```

### Full (8 columns):
```csv
source,target,score,category,distance,pval,size,type
cell1,cell2,0.85,ligand_a,12.5,0.001,15,typeA
cell2,cell3,0.72,ligand_b,8.3,0.003,22,typeB
```

**Just need:**
1. Source column (text)
2. Target column (text)
3. Weight column (number)
4. Optional: 5 more columns for visual features

**Column names don't matter** - you select them in the GUI!

---

## 🚀 Ready?

### Quick Start:
```bash
cd /app/examples
python load_network_examples.py
# Open MDV → Load "network_examples" → Visualize!
```

### Or Use Your Data:
```python
df = pd.read_csv("your_network.csv")
# Make sure it has: source, target, weight columns
# Load into MDV and visualize!
```

---

## 📚 More Help

- **Full Guide:** `/app/docs/FLEXIBLE_NETWORK_GUIDE.md`
- **Data Files:** `/app/examples/data/`
- **Python Examples:** `/app/python/mdvtools/examples/`

**Questions? Check the documentation or just try it - it's designed to be intuitive!**

---

**From CSV to interactive network in under 2 minutes!** 🎉🕸️

