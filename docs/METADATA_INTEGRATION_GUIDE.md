# Metadata Integration Guide

Guide for integrating external CSV/TSV metadata with MDV projects.

## Quick Decision Tree

```
Do you have the AnnData file?
├─ YES → Use Pre-Conversion Merge (Option 1)
│  └─ Easiest and cleanest approach
│
└─ NO → Is the MDV project already created?
   ├─ YES → Use CLI Tool (Option 2) or Python API (Option 3)
   └─ NO → Create project first, then add metadata
```

## Option 1: Pre-Conversion Merge (Recommended)

**When to use:** You have both the AnnData/Seurat file AND the metadata CSV before creating the MDV project.

**Advantages:**
- ✅ Cleanest approach
- ✅ All metadata included from the start
- ✅ Lineage tracks everything
- ✅ No need to modify MDV project after creation

### Example: Cell Phenotyping

```python
import scanpy as sc
import pandas as pd

# 1. Load your single-cell data
adata = sc.read_h5ad("experiment.h5ad")
print(f"Original: {adata.n_obs} cells")

# 2. Load external phenotype metadata
phenotypes = pd.read_csv("cell_phenotypes.csv")
# Assuming CSV has columns: cell_id, phenotype, confidence, annotator
print(f"Metadata: {len(phenotypes)} annotations")

# 3. Merge into adata.obs
adata.obs = adata.obs.merge(
    phenotypes,
    left_index=True,  # Use cell barcode as index
    right_on='cell_id',
    how='left'        # Keep all cells, even without phenotype
)

# 4. Save enriched data
adata.write_h5ad("experiment_with_phenotypes.h5ad")

# 5. Convert to MDV with all metadata
from mdvtools.umts.cli import convert_h5ad
```

Then convert:
```bash
python -m mdvtools.umts.cli convert-h5ad \
    experiment_with_phenotypes.h5ad \
    /app/mdv/experiment_project/
```

### Example: Multiple Metadata Sources

```python
import scanpy as sc
import pandas as pd

adata = sc.read_h5ad("data.h5ad")

# Load multiple metadata files
phenotypes = pd.read_csv("phenotypes.csv")
qc_metrics = pd.read_csv("qc_metrics.csv")
batch_info = pd.read_csv("batch_info.csv")

# Merge all metadata
metadata = phenotypes.merge(qc_metrics, on='cell_id', how='outer')
metadata = metadata.merge(batch_info, on='cell_id', how='outer')

# Add to AnnData
adata.obs = adata.obs.merge(
    metadata,
    left_index=True,
    right_on='cell_id',
    how='left'
)

adata.write_h5ad("data_enriched.h5ad")
```

## Option 2: CLI Tool (Post-Conversion)

**When to use:** MDV project already exists, need to add new metadata.

**Advantages:**
- ✅ Simple command-line interface
- ✅ No coding required
- ✅ Can be automated in scripts

### Usage

```bash
python -m mdvtools.umts.cli add-metadata \
    /app/mdv/my_project \
    cell_phenotypes.csv \
    --key cell_id \
    --datasource cells
```

### Full Example

```bash
# Your metadata CSV structure:
# cell_id,phenotype,confidence,method
# CELL_001,T-cell,0.95,manual
# CELL_002,B-cell,0.87,automated
# ...

# Add to existing project
python -m mdvtools.umts.cli add-metadata \
    /app/mdv/seurat-anndat-test \
    phenotypes.csv \
    --key cell_id \
    --datasource cells

# Result: Adds 'phenotype', 'confidence', 'method' columns to cells
```

### With Different Join Keys

```bash
# If your metadata uses 'barcode' instead of 'cell_id'
python -m mdvtools.umts.cli add-metadata \
    /app/mdv/my_project \
    metadata.tsv \
    --key cell_id \
    --metadata-key barcode \
    --datasource cells
```

### Overwrite Existing Columns

```bash
# Update existing annotations
python -m mdvtools.umts.cli add-metadata \
    /app/mdv/my_project \
    updated_phenotypes.csv \
    --key cell_id \
    --overwrite
```

## Option 3: Python API (Post-Conversion)

**When to use:** Need programmatic control, custom logic, or complex transformations.

### Basic Example

```python
from mdvtools import MDVProject
import pandas as pd

# Load project
mdv = MDVProject("/app/mdv/my_project")

# Load metadata
metadata = pd.read_csv("phenotypes.csv")

# Get current cell IDs
cell_ids = mdv.get_column("cells", "cell_id")

# Create lookup
phenotype_map = metadata.set_index('cell_id')['phenotype'].to_dict()

# Map values
phenotypes = [phenotype_map.get(cid, "Unknown") for cid in cell_ids]

# Add to project
mdv.set_column("cells", "phenotype", phenotypes)
```

### Advanced: Conditional Logic

```python
from mdvtools import MDVProject
import pandas as pd
import numpy as np

mdv = MDVProject("/app/mdv/my_project")

# Load metadata
metadata = pd.read_csv("phenotypes.csv")

# Get cell data
cell_ids = mdv.get_column("cells", "cell_id")
n_genes = mdv.get_column("cells", "n_genes")  # Assuming this exists

# Conditional phenotyping based on QC
phenotype_map = metadata.set_index('cell_id')['phenotype'].to_dict()
confidence_map = metadata.set_index('cell_id')['confidence'].to_dict()

phenotypes = []
for cid, ngenes in zip(cell_ids, n_genes):
    if ngenes < 200:
        # Low quality cells
        phenotypes.append("Low_quality")
    elif cid in phenotype_map and confidence_map.get(cid, 0) > 0.8:
        # High confidence annotation
        phenotypes.append(phenotype_map[cid])
    else:
        # Unknown or low confidence
        phenotypes.append("Unknown")

mdv.set_column("cells", "curated_phenotype", phenotypes)
```

### Batch Operations

```python
from mdvtools import MDVProject
import pandas as pd

def add_metadata_columns(project_path, metadata_file, join_key='cell_id'):
    """Add all columns from metadata file to project."""
    mdv = MDVProject(project_path)
    metadata = pd.read_csv(metadata_file)
    
    # Get join column
    cell_ids = mdv.get_column("cells", join_key)
    
    # Add each metadata column
    metadata_indexed = metadata.set_index(join_key)
    
    for col in metadata.columns:
        if col != join_key:
            values = [metadata_indexed.loc[cid, col] if cid in metadata_indexed.index 
                      else None for cid in cell_ids]
            mdv.set_column("cells", col, values)
            print(f"✓ Added {col}")

# Use it
add_metadata_columns("/app/mdv/my_project", "all_metadata.csv")
```

## When to Use MCP (Future Phase)

**MCP is NOT needed for simple metadata integration.** Use the options above.

**MCP would be useful for:**

### Scenario 1: Conflicting Annotations

```
Source A: CELL_001 → "T-cell" (confidence: 0.85)
Source B: CELL_001 → "NK-cell" (confidence: 0.92)

MCP Decision: Choose Source B (higher confidence)
or: Create "Ambiguous_T/NK" category
```

### Scenario 2: Schema Harmonization

```
File 1: cell_type → "CD4+ T cell"
File 2: phenotype → "T cell, CD4 positive"
File 3: annotation → "Tcell_CD4"

MCP: Recognizes these are the same thing → "CD4_T_cell"
```

### Scenario 3: Entity Resolution

```
File uses "Patient_001_Sample_A_CELL_0123"
MDV uses "CELL_0123"

MCP: Extracts barcode from complex ID
```

### Scenario 4: Multi-Modal Integration

```
scRNA-seq: 10,000 cells
ATAC-seq: 8,000 cells (subset)
Spatial: Different coordinate system

MCP: Aligns different modalities with different cell sets
```

## Common Use Cases

### Use Case 1: Flow Cytometry Gating

```bash
# You've done manual gating in FlowJo
# Export gate assignments as CSV

python -m mdvtools.umts.cli add-metadata \
    /app/mdv/facs_project \
    flowjo_gates.csv \
    --key cell_id
```

### Use Case 2: Automated Annotation Updates

```python
# Run automated annotation tool
import scanpy as sc
import celltypist

adata = sc.read_h5ad("data.h5ad")
predictions = celltypist.annotate(adata)

# Save predictions
predictions.to_csv("celltypist_annotations.csv")

# Add to MDV project
from mdvtools.umts.integrate_metadata import integrate_metadata
integrate_metadata(
    project_path="/app/mdv/my_project",
    metadata_file="celltypist_annotations.csv",
    join_key="cell_id"
)
```

### Use Case 3: Clinical Metadata

```bash
# Patient clinical data
# patient_id,age,sex,diagnosis,treatment

python -m mdvtools.umts.cli add-metadata \
    /app/mdv/patient_study \
    clinical_data.csv \
    --key patient_id \
    --datasource cells
```

### Use Case 4: Batch Effects

```csv
cell_id,batch,sequencing_date,library_prep,operator
CELL_001,batch1,2025-01-15,10x,Alice
CELL_002,batch1,2025-01-15,10x,Alice
CELL_003,batch2,2025-02-20,10x,Bob
```

```bash
python -m mdvtools.umts.cli add-metadata \
    /app/mdv/experiment \
    batch_info.csv
```

## Validation & QC

### Check Metadata Coverage

```python
from mdvtools import MDVProject
import pandas as pd

mdv = MDVProject("/app/mdv/my_project")
metadata = pd.read_csv("phenotypes.csv")

cell_ids = set(mdv.get_column("cells", "cell_id"))
metadata_ids = set(metadata['cell_id'])

print(f"Cells in MDV: {len(cell_ids)}")
print(f"Cells in metadata: {len(metadata_ids)}")
print(f"Matched: {len(cell_ids & metadata_ids)}")
print(f"Missing from metadata: {len(cell_ids - metadata_ids)}")
print(f"Extra in metadata: {len(metadata_ids - cell_ids)}")
```

### Verify Integration

```python
from mdvtools import MDVProject

mdv = MDVProject("/app/mdv/my_project")

# Check column was added
phenotypes = mdv.get_column("cells", "phenotype")
print(f"Phenotype column: {len(phenotypes)} values")
print(f"Unique phenotypes: {len(set(phenotypes))}")
print(f"Missing values: {phenotypes.count(None)}")

# Show distribution
from collections import Counter
print("\nPhenotype distribution:")
for pheno, count in Counter(phenotypes).most_common():
    print(f"  {pheno}: {count}")
```

## Summary

| Approach | When to Use | Complexity | Flexibility |
|----------|-------------|------------|-------------|
| **Pre-Conversion** | Have original data | ⭐ Simple | ⭐⭐⭐ High |
| **CLI Tool** | Post-conversion, simple | ⭐ Simple | ⭐⭐ Medium |
| **Python API** | Custom logic needed | ⭐⭐ Medium | ⭐⭐⭐ High |
| **MCP** | Complex conflicts | ⭐⭐⭐ Complex | ⭐⭐⭐ High |

**Recommendation:** Start with **Pre-Conversion** if possible, fall back to **CLI Tool** for post-hoc additions. Save MCP for future complex scenarios.

