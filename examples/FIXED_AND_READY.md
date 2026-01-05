# ✅ Network Examples - FIXED & READY

## ✨ Quick Demo (Working!)

```bash
cd /app/examples
python load_network_examples.py
```

**Output**:
```
✅ Loaded: ligand_receptor (30 rows, 8 columns)
✅ Loaded: simple (11 rows, 3 columns)  
✅ Loaded: genes (18 rows, 5 columns)
✅ PROJECT CREATED SUCCESSFULLY!
```

**Then**:
1. Open MDV GUI
2. Load project: `network_examples` (located in `/app/examples/`)
3. Visit views: "1. Ligand-Receptor", "2. Simple Network", "3. Gene Regulatory"
4. Follow instructions to add charts!

---

## 📁 CSV Files Available

All located in `/app/examples/data/`:

### 1. **ligand_receptor_network_example.csv** (30 rows)
```csv
source_cell,target_cell,interaction_score,ligand_type,spatial_distance,pvalue,cell_count,cell_type
T_cell_1,Endothelial_1,0.89,VEGFA,12.5,0.0001,18,T_cell
...
```

**All 8 columns** for complete visual encoding!

### 2. **simple_network_example.csv** (11 rows)
```csv
from_node,to_node,weight
A,B,0.85
...
```

**Minimal 3 columns** - easiest to understand!

### 3. **gene_regulatory_network_example.csv** (18 rows)
```csv
transcription_factor,target_gene,regulatory_score,regulation_type,binding_affinity
TP53,MDM2,0.95,activation,8.2
...
```

**Gene networks** with TF → gene interactions!

---

## 🎯 What Was Fixed

### Issue
- ❌ `add_datasource(data=df, size=len(df))` - Wrong parameters
- ❌ `project.save()` - Method doesn't exist

### Fixed
- ✅ `add_datasource(dataframe=df)` - Correct parameter
- ✅ Project saves automatically - no save() needed

### Files Updated
1. `/app/examples/load_network_examples.py` ✅
2. `/app/python/mdvtools/examples/flexible_network_example.py` ✅
3. `/app/python/mdvtools/examples/ligand_network_example.py` ✅

---

## 🚀 Quick Test

```bash
# Run the loader
cd /app/examples
python load_network_examples.py

# Should see:
# ✅ Loaded: ligand_receptor
# ✅ Loaded: simple
# ✅ Loaded: genes
# ✅ PROJECT CREATED SUCCESSFULLY!

# Project location: /app/examples/network_examples
```

---

## 📖 Using the Examples

### In MDV GUI:

1. **Load Project**:
   - File → Open Project
   - Navigate to `/app/examples/network_examples`
   - Or just select "network_examples" from recent

2. **Add Charts**:
   - View "2. Simple Network" (start here!)
   - Click "Add Chart"
   - Datasource: `simple`
   - Chart type: **"Network Graph (Flexible)"**
   - Select columns:
     * Source Node: `from_node`
     * Target Node: `to_node`
     * Edge Weight: `weight`
   - Add Chart → See network!

3. **Explore**:
   - Drag nodes
   - Click edges (highlights)
   - Settings (⚙️) → Adjust physics
   - Try other views!

---

## 💻 Use Your Own Data

```python
from mdvtools import MDVProject
import pandas as pd

# Load ANY CSV with network data
df = pd.read_csv("your_interactions.csv")

# Create project
project = MDVProject("my_network")

# Add datasource - CORRECT SYNTAX:
project.add_datasource("data", dataframe=df)
# NOT: data=df, size=len(df) ❌

# No need to call save() - automatic!
```

Then open in GUI and visualize!

---

## ✅ Status

- Scripts: **WORKING** ✅
- CSV files: **READY** ✅  
- Documentation: **COMPLETE** ✅
- Examples: **TESTED** ✅

**Everything is ready to use!** 🎉

---

## 📚 Full Documentation

- **Quick Start**: `/app/examples/NETWORK_EXAMPLES_QUICKSTART.md`
- **CSV Info**: `/app/examples/data/README.md`
- **Chart Guide**: `/app/docs/FLEXIBLE_NETWORK_GUIDE.md`
- **Full Summary**: `/app/docs/NETWORK_CHARTS_FINAL_SUMMARY.md`

---

**Ready to visualize networks in MDV!** 🕸️📊

