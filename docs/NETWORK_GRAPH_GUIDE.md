# Network Graph - Interactive Network Visualization

The **Network Graph** creates interactive force-directed network visualizations that work with **any table**. No Python configuration required - just select your columns in the GUI and visualize!

**Chart Type**: `flexible_network_chart`  
**Display Name**: "Network Graph"  
**Location**: `/app/src/charts/FlexibleNetworkChart.js`

---

## Quick Start (2 Minutes)

### Your Data (any CSV/table):
```csv
from_node,to_node,weight,category,distance,pvalue,count,type
A,B,0.85,ligand,12.5,0.001,15,T_cell
B,C,0.72,ligand,8.3,0.003,22,Cancer
C,D,0.91,receptor,15.2,0.0001,8,Fibroblast
```

### Load Data (Python):
```python
from mdvtools import MDVProject
import pandas as pd

df = pd.read_csv("network.csv")
project = MDVProject("my_network")
project.add_datasource("data", dataframe=df)
project.save()
```

### Add Chart (GUI):
1. Open MDV → Load project
2. Click **"Add Chart"**
3. Select datasource: `data`
4. Chart type: **"Network Graph"**
5. **Select columns**:
   - Source Node → `from_node` ✅ Required
   - Target Node → `to_node` ✅ Required
   - Edge Weight → `weight` ✅ Required
   - Category Filter → `category` ⭕ Optional
   - Edge Length → `distance` ⭕ Optional
   - Edge Color → `pvalue` ⭕ Optional
   - Node Size → `count` ⭕ Optional
   - Node Type → `type` ⭕ Optional
6. Click **"Add Chart"**

**Done! Interactive network appears immediately!** 🎉

---

## Column Guide

### Required Columns (3)

| Column | Type | Purpose | Example Values |
|--------|------|---------|----------------|
| **Source Node** | Text | Start of connection | `cell_1`, `gene_A`, `user_alice` |
| **Target Node** | Text | End of connection | `cell_2`, `gene_B`, `user_bob` |
| **Edge Weight** | Number | Strength/score | `0.85`, `42`, `0.95` |

### Optional Columns (5)

| Column | Type | Visual Property | Example Values |
|--------|------|----------------|----------------|
| **Category Filter** | Text | Filter dropdown | `VEGFA`, `pathway_1`, `ligand` |
| **Edge Length** | Number | Link distance in layout | `12.5`, `8.3`, `100` |
| **Edge Color** | Number | Link color gradient | `0.001`, `0.05`, `1.0` |
| **Node Size** | Number | Node circle radius | `1`, `50`, `1000` |
| **Node Type** | Text | Node color categories | `T_cell`, `Cancer`, `admin` |

---

## Visual Encoding

| Visual Property | Data Mapping | Example |
|----------------|--------------|---------|
| **Link Thickness** | Edge Weight column | 0.1 → thin line, 1.0 → thick line |
| **Link Length** | Edge Length column | Small value → nodes closer, Large value → nodes farther |
| **Link Color** | Edge Color column | Low → red, Medium → orange, High → gray |
| **Node Radius** | Node Size column | Small value → small circle, Large value → big circle |
| **Node Color** | Node Type column | Each category gets unique color |

---

## Interactive Features

### Mouse Controls
- **Scroll wheel/trackpad** → Zoom in/out
- **Click + drag** → Pan the network
- **Drag node** → Reposition manually (stays pinned)
- **Click link** → Highlight that connection
- **Click node** → Highlight all connected links
- **Hover link** → Show tooltip with details

### Settings Panel (⚙️)

#### Physics Controls
- **Link Strength** (0-2): How strongly links pull nodes together
  - Low (0.3) = Loose, spread out network
  - High (1.5) = Tight, compact network
- **Node Repulsion** (-1000 to 0): How strongly nodes push apart
  - Strong (-800) = Very spread out
  - Weak (-100) = Nodes can be close

#### Display Controls
- **Node Radius** (3-30): Base size for all nodes (if no Node Size column)
- **Show Node Labels** (✓/✗): Display node IDs as text
- **Label Size** (5-20): Font size for labels
- **Link Opacity** (0.1-1.0): Transparency of links
- **Show Directionality** (✓/✗): Display arrows on links

#### Category Filter (if Category column selected)
- **Dropdown**: Select which category to display
- **All**: Show all categories

#### Zoom
- **Reset Zoom** button: Return to default view (animated)

#### Color Settings (if "Color By" selected)
- **Color By**: Select numeric column for node coloring
- **Show Color Legend**: Toggle color scale display
- **SymLog Color Scale**: Use logarithmic color mapping
- **Trim Color Scale**: Remove outliers from color range
- **Treat zero as missing**: Handle zero values specially

#### Parameters (Advanced)
- Change any source/target/weight/optional columns
- Network rebuilds automatically

---

## Common Use Cases

### 1. Ligand-Receptor Interactions
```csv
ligand,source_cell,target_cell,score,pvalue,distance,count,cell_type
VEGFA,T_cell_1,Endothelial_1,0.85,0.001,12.5,15,T_cell
TGFB1,Fibroblast_1,Cancer_1,0.91,0.0001,15.2,8,Fibroblast
```
→ Visualize cell-cell communication

### 2. Gene Regulatory Networks
```csv
tf,gene,regulation_score
TP53,MDM2,0.95
STAT3,IL6,0.92
NF-KB,TNF,0.88
```
→ Understand transcription factor networks

### 3. Social Networks
```csv
user_a,user_b,connection_strength
alice,bob,0.95
bob,charlie,0.82
charlie,david,0.77
```
→ Analyze user connections

### 4. Protein-Protein Interactions
```csv
protein_a,protein_b,confidence,experiment,binding_affinity
BRCA1,BARD1,0.98,Y2H,50
TP53,MDM2,0.95,CoIP,75
```
→ Study molecular interactions

---

## Tips & Tricks

### Making Networks More Readable

**Spread out dense networks**:
- Link Strength: 0.3
- Node Repulsion: -800
- Link Opacity: 0.4

**Compact sparse networks**:
- Link Strength: 1.5
- Node Repulsion: -100

**Focus on important connections**:
1. Use Edge Weight to size by importance
2. Use Edge Color to highlight significance
3. Lower Link Opacity to fade weak connections

**Find specific nodes**:
1. Use Category Filter to show subset
2. Click node to highlight its connections
3. Use Node Type colors to identify groups

### Handling Large Networks

**For 100+ nodes**:
- Increase Node Repulsion (-800 to -1000)
- Reduce Link Opacity (0.2-0.4)
- Turn off labels initially
- Use Category Filter to show subsets
- Zoom in on areas of interest

**For dense connections**:
- Use Edge Color to show significance
- Filter to top N% of edges by weight
- Use directionality to show flow

### Performance

- Chart handles 1000+ nodes smoothly
- Force simulation stabilizes in 1-2 seconds
- Zoom/pan is instant
- Dragging individual nodes is responsive

---

## Example Workflows

### Explore CellPhoneDB Results
```python
# Load CellPhoneDB output
means = pd.read_csv("means.txt", sep="\t")

# Transform to network format
df = means.melt(id_vars=['ligand', 'receptor'])
df.columns = ['ligand', 'receptor', 'cell_pair', 'score']
df[['source', 'target']] = df['cell_pair'].str.split('|', expand=True)

# Load in MDV
project = MDVProject("cellphonedb_viz")
project.add_datasource("interactions", dataframe=df)
project.save()
```
→ Select columns in GUI, filter by ligand, explore!

### Compare Network States
```python
# Load two conditions
df_before = pd.read_csv("network_before.csv")
df_after = pd.read_csv("network_after.csv")

project = MDVProject("comparison")
project.add_datasource("before", dataframe=df_before)
project.add_datasource("after", dataframe=df_after)
project.save()
```
→ Add both as separate charts, compare visually

### Integrate with Spatial Data
```python
# Calculate distances between cells
interactions['distance'] = calculate_spatial_distance(
    interactions['source_cell'],
    interactions['target_cell'],
    spatial_coords
)

project.add_datasource("spatial_network", dataframe=interactions)
```
→ Use distance as Edge Length, shows spatial relationships!

---

## Troubleshooting

| Problem | Solution |
|---------|----------|
| Network too spread out | Increase Link Strength (try 1.0-1.5) |
| Network too compact | Increase Node Repulsion (try -600 to -800) |
| Can't see labels | Increase Label Size or zoom in |
| Links too faint | Increase Link Opacity |
| Too many nodes | Use Category Filter to show subset |
| Wrong connections showing | Check Source/Target columns selected correctly |
| Self-loops visible | Chart automatically skips self-loops (A→A) |
| Bidirectional links | Chart keeps strongest direction only |
| Colors all wrong | Check Node Type column has correct categories |
| Zoom not working | Refresh page - zoom should work immediately |

---

## Advanced Features

### Bidirectional Link Handling
When both A→B and B→A exist, the chart:
1. Detects both directions
2. Compares Edge Weight values
3. Keeps only the stronger direction
4. Shows it as single link

### Filtered Node Indication
Nodes with hidden connections (due to MDV filters):
- Display at 40% opacity
- All visible connections shown normally
- Indicates incomplete network view

### Color By Any Numeric Column
Beyond the 8 basic columns:
1. Open Settings → Color Settings
2. Select "Color By" → Choose ANY numeric column
3. Adjust color scale, trim outliers, use log scale
4. Color legend shows automatically

### Parameter Changes
Change any column mapping:
1. Settings → Parameters
2. Select different column from dropdown
3. Network rebuilds automatically
4. Maintains zoom/pan state

---

## Example CSV Files

See `/app/examples/data/` for ready-to-use examples:

1. **`ligand_receptor_network_example.csv`** - Full 8-column demo
2. **`simple_network_example.csv`** - Minimal 3-column demo  
3. **`gene_regulatory_network_example.csv`** - TF network demo

Load with:
```bash
python /app/examples/load_network_examples.py
```

---

## Comparison with Other Network Tools

### vs. Cytoscape
- ✅ Faster setup (no export/import)
- ✅ Integrated with MDV data
- ✅ Real-time filtering
- ⭕ Fewer layout algorithms (only force-directed)

### vs. NetworkX/igraph
- ✅ Interactive GUI (not static images)
- ✅ No coding required
- ✅ Live updates
- ⭕ Fewer graph metrics

### vs. Gephi
- ✅ No separate application
- ✅ Web-based
- ✅ Part of analysis workflow
- ⭕ Fewer styling options

**Best for**: Exploratory analysis, quick visualization, integrated with other MDV charts

---

## Technical Details

- **D3.js force simulation**: Fast, smooth physics
- **SVG rendering**: Crisp visuals at any zoom
- **D3 zoom behavior**: Native pan/zoom support
- **Automatic scale generation**: Quantile-based for robustness
- **Filter integration**: Respects MDV global filters
- **Data highlighting**: Bidirectional selection with other charts

---

## Files

- **Chart**: `/app/src/charts/FlexibleNetworkChart.js`
- **Example Script**: `/app/examples/load_network_examples.py`
- **Example Data**: `/app/examples/data/*network*.csv`
- **This Guide**: `/app/docs/NETWORK_GRAPH_GUIDE.md`

---

## Quick Reference

**Minimum viable network**: 3 columns (source, target, weight)  
**Maximum features**: 8 columns (all visual properties)  
**Setup time**: 2 minutes  
**Python config**: None required  
**Chart updates**: Automatic on filter/parameter change  
**Mouse**: Scroll to zoom, drag to pan, click to interact  
**Reset**: Settings → Reset Zoom button  

**Works with ANY table format - just select your columns!** 🚀


