# Network Visualization in MDV - Complete Solution ✅

## Overview

MDV now has **TWO network chart options** to fit different workflows:

| Chart | Best For | Setup | Friction |
|-------|----------|-------|----------|
| **Network Graph (Flexible)** 🆕 | Exploration, varied data | None - GUI only | ⭐ Zero |
| **Spatial Connectivity Map** | Production, standardized data | Python metadata | Medium |

**Both are now available in the GUI!**

---

## Option 1: Network Graph (Flexible) 🆕 **[RECOMMENDED FOR MOST USERS]**

### ✨ The Zero-Friction Solution

**Perfect for users who want:**
- ✅ No Python configuration
- ✅ Works with ANY table format
- ✅ Direct column selection in GUI
- ✅ Quick exploration
- ✅ Maximum flexibility

### Quick Start (2 minutes)

**Your data** (any format):
```csv
from,to,weight,category,distance,pvalue,size,type
A,B,0.85,ligand,12.5,0.001,15,T_cell
B,C,0.72,ligand,8.3,0.003,22,Cancer
```

**Python** (just to load):
```python
from mdvtools import MDVProject
import pandas as pd

df = pd.read_csv("network.csv")
project = MDVProject("my_network")
project.add_datasource("data", data=df, size=len(df))
project.save()
```

**GUI**:
1. Open MDV, load project
2. Click **"Add Chart"**
3. Select datasource: `data`
4. Chart type: **"Network Graph (Flexible)"**
5. **Select columns**:
   - Source Node → `from` ✅ Required
   - Target Node → `to` ✅ Required
   - Edge Weight → `weight` ✅ Required
   - Category Filter → `category` ⭕ Optional
   - Edge Length → `distance` ⭕ Optional
   - Edge Color → `pvalue` ⭕ Optional
   - Node Size → `size` ⭕ Optional
   - Node Type → `type` ⭕ Optional
6. Click **"Add Chart"** → Done! 🎉

**Result**: Interactive network with all visual encodings!

### Files
- **Chart**: `/app/src/charts/FlexibleNetworkChart.js`
- **Guide**: `/app/docs/FLEXIBLE_NETWORK_GUIDE.md`
- **Example**: `/app/python/mdvtools/examples/flexible_network_example.py`

---

## Option 2: Spatial Connectivity Map (Existing)

### 🏗️ The Production-Ready Solution

**Perfect for users who want:**
- ✅ Validated, consistent configurations
- ✅ Standardized data pipelines
- ✅ Metadata-driven approach
- ✅ Integration with spatial tools

### Quick Start (Python setup required)

**Python** (configuration):
```python
from mdvtools import MDVProject, setup_ligand_network
import pandas as pd

df = pd.read_csv("interactions.csv")
project = MDVProject("configured_network")
project.add_datasource("interactions", data=df, size=len(df))

# Configure metadata
setup_ligand_network(
    project,
    "interactions",
    ligand_column="ligand_type",
    source_cell_column="source_cell",
    target_cell_column="target_cell",
    interaction_score_column="score",
    # ... more options
)

project.save()
```

**GUI**:
1. Open MDV, load project
2. Click **"Add Chart"**
3. Select datasource: `interactions`
4. Chart type: **"Spatial Connectivity Map"** ← Only appears if metadata configured
5. **Select category**: Choose ligand type from dropdown
6. Click **"Add Chart"** → Done!

**Result**: Pre-configured network with validated columns!

### Files
- **Chart**: `/app/src/charts/CellNetworkChart.js`
- **Helper**: `/app/python/mdvtools/network_helpers.py`
- **Guide**: `/app/docs/NETWORK_CHART_GUIDE.md`
- **Example**: `/app/python/mdvtools/examples/ligand_network_example.py`

---

## Feature Comparison

| Feature | Flexible Network | Spatial Connectivity |
|---------|-----------------|---------------------|
| **Python setup** | None | Required |
| **Metadata required** | No | Yes (`interactions`) |
| **Works with** | Any table | Configured datasources only |
| **Column selection** | 8 GUI dropdowns | Pre-configured |
| **Minimum columns** | 3 | 4+ (via metadata) |
| **GUI clicks to add** | 6-8 | 3-4 |
| **Learning curve** | Low | Medium |
| **Flexibility** | Maximum | Limited |
| **Validation** | None | Column type checking |
| **Use case** | Exploration | Production |
| **Physics controls** | ✅ Yes | ✅ Yes |
| **Visual encoding** | ✅ 5 channels | ✅ 5 channels |
| **Bidirectional handling** | ✅ Yes | ✅ Yes |
| **Ego network mode** | ❌ No | ✅ Yes |
| **Spatial linking** | ⭕ Manual | ✅ Built-in |

---

## When to Use Each

### Use **Flexible Network** When:

✅ **You want to explore quickly**
- No time for Python configuration
- Just want to see the data as a network
- Trying different column combinations

✅ **Data format varies**
- Different tables have different column names
- Working with many different datasets
- Non-standard data structures

✅ **You're not a Python expert**
- Prefer GUI-only workflows
- Don't want to write code
- Want immediate results

✅ **Rapid prototyping**
- Testing ideas quickly
- Showing quick visualizations
- Iterating on analysis

### Use **Spatial Connectivity Map** When:

✅ **You have standardized pipelines**
- Same data format every time
- Production workflows
- Automated analysis

✅ **You want validation**
- Ensure columns are correct types
- Prevent configuration errors
- Team consistency

✅ **You need advanced features**
- Ego network mode (focus on one cell)
- Complex spatial integration
- Multi-level network exploration

✅ **Building production tools**
- Delivering to non-technical users
- Want pre-configured, validated setup
- Part of larger workflow

---

## Visual Encoding (Both Charts)

| Property | Maps To | Example |
|----------|---------|---------|
| **Link Thickness** | Edge strength/score | 0.1 → thin, 1.0 → thick |
| **Link Length** | Spatial distance or score | 5µm → short, 50µm → long |
| **Link Color** | Significance (p-value) | 0.001 → red, 0.05 → orange, 1.0 → gray |
| **Node Size** | Count/importance | 1 → small, 100 → large |
| **Node Color** | Type/category | T_cell → blue, Cancer → red |

---

## Common Use Cases

### 1. Ligand-Receptor Interactions

**Flexible** (Quick exploration):
```python
# Just load CSV
df = pd.read_csv("cellphonedb_means.csv")
project.add_datasource("lr", data=df, size=len(df))
# Select columns in GUI!
```

**Spatial Connectivity** (Production):
```python
# Configure once
setup_ligand_network(project, "lr", 
                     ligand_column="ligand",
                     source_cell_column="sender",
                     target_cell_column="receiver")
# Users just select ligand type in GUI
```

### 2. Gene Regulatory Networks

**Flexible**:
```csv
tf,gene,score
TP53,MDM2,0.95
STAT3,IL6,0.92
```
→ Load and select in GUI

### 3. Cell Communication

**Flexible**:
```csv
from_cell,to_cell,pathway,score
Neuron_1,Astrocyte_1,glutamate,0.88
```
→ Load and select in GUI

### 4. Social Networks

**Flexible** (works with ANY graph!):
```csv
user_a,user_b,connection_strength
alice,bob,0.95
bob,charlie,0.82
```

---

## Interactive Features (Both Charts)

### Mouse Controls
- **Drag nodes**: Manually reposition
- **Click link**: Highlight edge
- **Click node**: Highlight all connected edges
- **Hover**: Tooltips with details

### Physics Controls (⚙️ Settings)
- **Link Strength** (0-2): Tightness of connections
- **Node Repulsion** (-1000 to 0): Spacing between nodes

### Display Controls
- **Node Labels**: Show/hide text
- **Directionality**: Show/hide arrows
- **Node Size**: Adjust radius
- **Color Schemes**: Choose visual style

---

## Integration Examples

### From CellPhoneDB
```python
# Works with either chart!
means = pd.read_csv("means.txt", sep="\t")
# Reshape to network format...
project.add_datasource("cellphonedb", data=df, size=len(df))

# Flexible: Select columns in GUI
# Spatial: setup_ligand_network()
```

### From NicheNet / LIANA / CellChat
```python
# Any cell communication tool output works!
df = load_results("tool_output.csv")
project.add_datasource("comm", data=df, size=len(df))
# Choose chart type based on workflow
```

---

## Files Created/Modified

### New Files 🆕
1. `/app/src/charts/FlexibleNetworkChart.js` - New flexible chart
2. `/app/python/mdvtools/network_helpers.py` - Helper functions
3. `/app/docs/FLEXIBLE_NETWORK_GUIDE.md` - Flexible chart guide
4. `/app/docs/NETWORK_CHART_GUIDE.md` - Spatial connectivity guide
5. `/app/docs/NETWORK_CHART_CHEATSHEET.md` - Quick reference
6. `/app/python/mdvtools/examples/flexible_network_example.py` - Flexible example
7. `/app/python/mdvtools/examples/ligand_network_example.py` - Configured example

### Modified Files ✏️
1. `/app/src/charts/ChartManager.js` - Added FlexibleNetworkChart import (line 41)
2. `/app/python/mdvtools/__init__.py` - Added network helper exports

### Existing Files (No Changes) ✅
1. `/app/src/charts/CellNetworkChart.js` - Already working!

---

## Testing Checklist

### Flexible Network Chart
- [x] JavaScript syntax valid
- [x] Registered in ChartManager
- [x] No metadata required
- [x] 8 parameter dropdowns defined
- [x] Example created
- [x] Documentation complete

### Spatial Connectivity Map
- [x] Already working (no changes)
- [x] Python helpers created
- [x] Example updated
- [x] Documentation added

### Integration
- [x] Both charts in ChartManager
- [x] Different chart type names
- [x] Module exports working
- [x] Examples compile successfully

---

## User Decision Tree

```
Do you have interaction/network data?
  ↓
  YES
  ↓
Do you want to configure in Python or just use the GUI?
  ↓
  ┌─────────┬─────────┐
  GUI ONLY  │  PYTHON │
  ↓         │    ↓    │
  Use       │   Use   │
  Flexible  │ Spatial │
  Network   │ Connect.│
```

Or simply: **Try Flexible first - it's easier!**

---

## Quick Commands

### Run Flexible Example:
```bash
python /app/python/mdvtools/examples/flexible_network_example.py
```

### Run Configured Example:
```bash
python /app/python/mdvtools/examples/ligand_network_example.py
```

### Import in Python:
```python
# Flexible - just load data
from mdvtools import MDVProject

# Configured - use helpers
from mdvtools import setup_ligand_network
```

---

## Summary

**Mission Accomplished! ✅**

You now have **TWO powerful network visualization options**:

1. **🆕 Network Graph (Flexible)** - Zero friction, maximum flexibility
   - No setup, works with any table
   - Perfect for exploration

2. **🏗️ Spatial Connectivity Map** - Production-ready, validated
   - Python configuration, consistent results
   - Perfect for pipelines

**Both charts offer:**
- Interactive force-directed layouts
- Rich visual encoding (5 channels)
- Physics controls
- Integration with MDV ecosystem
- Bidirectional link handling

**Choose based on your workflow - both are excellent!** 🎉

---

## Next Steps

1. **Try the Flexible chart first** (easiest)
   - Run: `python flexible_network_example.py`
   - Open in GUI
   - Add chart with column selection

2. **Explore Spatial Connectivity** (if needed)
   - Run: `python ligand_network_example.py`
   - See metadata configuration
   - Compare workflows

3. **Use with your data**
   - Load your network data
   - Choose chart type
   - Visualize!

**Happy networking!** 🕸️📊🚀

