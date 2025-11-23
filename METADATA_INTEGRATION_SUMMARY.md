# Metadata Integration - Quick Summary

## ✅ You Now Have 3 Ways to Add External Metadata

### 1. Pre-Conversion (Best for New Projects)

```python
import scanpy as sc
import pandas as pd

# Load data and metadata
adata = sc.read_h5ad("data.h5ad")
metadata = pd.read_csv("phenotypes.csv")

# Merge before conversion
adata.obs = adata.obs.merge(metadata, left_index=True, right_on='cell_id', how='left')
adata.write_h5ad("enriched.h5ad")

# Convert to MDV
```

```bash
python -m mdvtools.umts.cli convert-h5ad enriched.h5ad /app/mdv/project/
```

### 2. CLI Tool (Best for Existing Projects)

```bash
# Add cell phenotypes from CSV
python -m mdvtools.umts.cli add-metadata \
    /app/mdv/seurat-anndat-test \
    cell_phenotypes.csv \
    --key cell_id

# Add with different column names
python -m mdvtools.umts.cli add-metadata \
    /app/mdv/my_project \
    metadata.tsv \
    --key cell_id \
    --metadata-key barcode

# Update existing columns
python -m mdvtools.umts.cli add-metadata \
    /app/mdv/my_project \
    updated_data.csv \
    --overwrite
```

### 3. Python API (Best for Custom Logic)

```python
from mdvtools import MDVProject
import pandas as pd

mdv = MDVProject("/app/mdv/my_project")
metadata = pd.read_csv("phenotypes.csv")

# Get cells
cell_ids = mdv.get_column("cells", "cell_id")

# Create mapping
phenotype_map = metadata.set_index('cell_id')['phenotype'].to_dict()
phenotypes = [phenotype_map.get(cid, "Unknown") for cid in cell_ids]

# Add to project
mdv.set_column("cells", "phenotype", phenotypes)
```

## When to Use Each Approach

| Situation | Use This | Why |
|-----------|----------|-----|
| Converting new data | Pre-Conversion | Cleanest, all tracked in lineage |
| Existing MDV project | CLI Tool | Fast, no coding |
| Complex logic needed | Python API | Full flexibility |
| Conflicting sources | Future MCP | Intelligent resolution |

## MCP is NOT Needed For Simple Cases

**You DON'T need MCP for:**
- ✅ Adding CSV/TSV metadata to cells
- ✅ Merging phenotype annotations
- ✅ Adding batch information
- ✅ Integrating QC metrics

**MCP would help with (Future):**
- ❌ Conflicting annotations from multiple sources
- ❌ Different ID formats that need reconciliation
- ❌ Schema harmonization across datasets
- ❌ Multi-modal alignment with different cell sets

## Your Workflow

For your current use case (adding phenotyping metadata):

```bash
# Option A: Before conversion
# 1. Merge in Python/R
# 2. Convert to MDV

# Option B: After conversion (RECOMMENDED FOR YOU)
python -m mdvtools.umts.cli add-metadata \
    /app/mdv/seurat-anndat-test \
    your_phenotypes.csv \
    --key cell_id
```

Then refresh the web interface:
```
http://localhost:5055/rescan_projects
```

## Files Created

- `/app/python/mdvtools/umts/integrate_metadata.py` - Core integration logic
- `/app/python/mdvtools/umts/cli.py` - Updated with `add-metadata` command
- `/app/docs/METADATA_INTEGRATION_GUIDE.md` - Full documentation

## Quick Test

Create a test metadata file:

```bash
cat > /tmp/test_metadata.csv << EOF
cell_id,phenotype,confidence
CELL_001,T-cell,0.95
CELL_002,B-cell,0.87
CELL_003,Monocyte,0.92
EOF

# Add to your project
python -m mdvtools.umts.cli add-metadata \
    /app/mdv/seurat-anndat-test \
    /tmp/test_metadata.csv \
    --key cell_id
```

## Summary

✅ **Simple metadata integration: Use CLI or Python API**  
❌ **MCP: Only for complex future scenarios**  
📝 **Documentation: See `/app/docs/METADATA_INTEGRATION_GUIDE.md`**  

Your phenotyping CSV can be integrated right now with a single command!

