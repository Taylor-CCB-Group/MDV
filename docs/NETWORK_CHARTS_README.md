# Network Visualization in MDV

## Network Graph - Interactive Network Visualization

**Zero-friction network visualization that works with any table** ⭐

- ✅ **No Python configuration required**
- ✅ **Select columns directly in GUI** (8 intuitive dropdowns)
- ✅ **Works with ANY table format** (no metadata needed)
- ✅ **Zoom, pan, color by column, physics controls**
- ✅ **Perfect for exploration and analysis**

**[→ Full Guide: NETWORK_GRAPH_GUIDE.md](NETWORK_GRAPH_GUIDE.md)**

### Quick Start (2 minutes)
```python
# Just load your data
import pandas as pd
from mdvtools import MDVProject

df = pd.read_csv("network.csv")  # Any CSV with from/to/weight columns
project = MDVProject("my_network")
project.add_datasource("data", dataframe=df)
project.save()
```

**Then in GUI:**
1. Click "Add Chart"
2. Select datasource
3. Chart type: **"Network Graph"**
4. Select columns from dropdowns:
   - Source Node (required)
   - Target Node (required)
   - Edge Weight (required)
   - 5 optional columns for visual properties
5. Done! Interactive network appears

**That's it - no configuration, no metadata, just works!** 🚀

---

## Advanced: Spatial Connectivity Map

For users who need advanced features like ego networks and spatial integration, MDV also includes the **Spatial Connectivity Map** chart. This requires Python metadata configuration and is only visible in the GUI after setup.

**Most users should use Network Graph** - it's simpler and more flexible. Only use Spatial Connectivity Map if you specifically need its advanced features.

[→ Learn about Spatial Connectivity Map](NETWORK_CHART_GUIDE.md) (Advanced users only)

---

## Features

- **Interactive force-directed layout** with physics controls
- **Zoom and pan** with mouse/trackpad
- **8 visual encoding channels**: link thickness, length, color, node size, color, and more
- **Color by any numeric column** with automatic legends
- **Category filtering** with GUI dropdown
- **Settings panel** for physics, display, and parameter changes
- **Drag nodes** to reposition manually
- **Click to highlight** connections
- **Bidirectional link handling** (keeps strongest direction)
- **Filter integration** with MDV ecosystem

---

## Example Files

- **Example Data**: `/app/examples/data/*network*.csv`
- **Quick Load**: `python /app/examples/load_network_examples.py`

---

## Documentation

**→ [Network Graph Complete Guide](NETWORK_GRAPH_GUIDE.md)** ⭐ **(Start here!)**

[Spatial Connectivity Map Guide](NETWORK_CHART_GUIDE.md) (Advanced - requires Python setup)


