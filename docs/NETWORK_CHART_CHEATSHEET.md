# Network Chart Quick Reference 🚀

## 30-Second Setup

```python
from mdvtools import MDVProject, setup_ligand_network
import pandas as pd

df = pd.read_csv("interactions.csv")
project = MDVProject("network_viz")
project.add_datasource("interactions", data=df, size=len(df))

setup_ligand_network(project, "interactions", 
                     ligand_column="ligand_type",
                     source_cell_column="source_cell", 
                     target_cell_column="target_cell",
                     interaction_score_column="score")
project.save()
```

Then in GUI: **Add Chart → "Spatial Connectivity Map" → Select ligand → Done!**

---

## Required Data Columns

| Column | Type | Example |
|--------|------|---------|
| Category/Ligand | Text | "VEGFA", "TGFB1" |
| Source Cell | Text | "T_cell_1" |
| Target Cell | Text | "Endothelial_1" |
| Score | Number | 0.85 |

---

## Visual Encoding

| Property | Shows |
|----------|-------|
| **Link thickness** | Interaction strength |
| **Link color** | P-value/significance |
| **Node size** | Cell count |
| **Node color** | Cell type |

---

## Common Settings

```javascript
// Spread out network
Link Strength: 0.3
Node Repulsion: -800

// Compact network  
Link Strength: 1.5
Node Repulsion: -100

// Ego network (focus on one cell)
Center Cell: "T_cell_1"
Levels: 2
```

---

## Keyboard Shortcuts

| Action | How |
|--------|-----|
| Drag node | Click + drag |
| Highlight | Click link/node |
| Settings | ⚙️ icon |
| Zoom | Mouse wheel |

---

## Integration Examples

### CellPhoneDB
```python
df = pd.read_csv("means.txt", sep="\t")
# Transform to: ligand_type, source_cell_id, target_cell_id, score
setup_ligand_network(project, "cellphonedb", ...)
```

### Spatial Link
```python
setup_ligand_network(project, "interactions",
                     cells_datasource="cells",  # ← Links to spatial
                     spatial_distance_column="distance")
```

---

## Troubleshooting

| Problem | Solution |
|---------|----------|
| Chart doesn't appear | Run `setup_ligand_network()` |
| Empty network | Select ligand type in dropdown |
| Too spread out | Increase link strength |
| Too compact | Decrease node repulsion |

---

## Files

- **Guide**: `/app/docs/NETWORK_CHART_GUIDE.md`
- **Example**: `/app/python/mdvtools/examples/ligand_network_example.py`
- **API**: `/app/python/mdvtools/network_helpers.py`

---

**That's it! Network in 30 seconds. Full guide available above.**

