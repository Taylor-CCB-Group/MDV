# UMTS Command Line Interface

The UMTS CLI provides command-line tools for converting various single-cell data formats to MDV projects with automatic lineage tracking.

## Installation

### Basic Installation

```bash
pip install mdvtools
```

This provides lineage tracking for h5ad files.

### Full Installation (with Seurat Support)

To enable Seurat RDS file conversion, install the optional dependencies:

```bash
# Install Python packages
pip install mdvtools[seurat]

# Install R packages (in R console)
install.packages("Seurat")
remotes::install_github("mojaveazure/seurat-disk")
```

**Note:** Seurat conversion requires R to be installed on your system.

## Commands

### Check Dependencies

Verify that all required dependencies are installed:

```bash
umts check-deps
```

**Output:**
```
Checking UMTS dependencies...

Python Packages:
  ✓ rpy2
  ✓ anndata2ri

R Environment:
  ✓ R available
  ✓ Seurat package
  ✓ SeuratDisk package

Overall Status: ✓ Ready
```

### Convert Seurat RDS to MDV

Convert a Seurat object from an RDS file to MDV format:

```bash
umts convert-seurat input.rds output_project/
```

**Options:**
- `--max-dims N` - Maximum dimensions to include (default: 3)
- `--label PREFIX` - Prefix for datasource names
- `--chunk-data` - Process data in chunks (slower, lower memory)
- `--no-delete` - Merge with existing project instead of overwriting
- `--no-lineage` - Disable lineage tracking
- `-v, --verbose` - Show detailed error messages

**Example:**
```bash
# Basic conversion
umts convert-seurat pbmc_seurat.rds ./pbmc_mdv/

# With custom options
umts convert-seurat pbmc_seurat.rds ./pbmc_mdv/ \
    --max-dims 2 \
    --label rna_ \
    --chunk-data
```

**Output:**
```
Checking dependencies...
✓ All dependencies available

Loading Seurat object from pbmc_seurat.rds...
Converting Seurat to h5ad format...
Loading as AnnData...
Loaded: 2638 cells, 1838 genes
Converting to MDV project at ./pbmc_mdv/...
Getting Matrix
Adding gene expression

✓ Successfully created MDV project at: ./pbmc_mdv/
✓ Lineage information saved
  - Source: /path/to/pbmc_seurat.rds
  - SHA256: a1b2c3d4e5f6...
  - Timestamp: 2025-11-23T15:30:00Z
```

### Convert h5ad to MDV

Convert an AnnData h5ad file to MDV format:

```bash
umts convert-h5ad input.h5ad output_project/
```

**Options:**
Same as `convert-seurat` command.

**Example:**
```bash
# Basic conversion
umts convert-h5ad pbmc3k.h5ad ./pbmc3k_mdv/

# With custom options
umts convert-h5ad pbmc3k.h5ad ./pbmc3k_mdv/ \
    --max-dims 3 \
    --no-lineage
```

### Show Lineage Information

Display lineage and provenance information for an existing MDV project:

```bash
umts show-lineage project_directory/
```

**Options:**
- `--packages` - Show all package versions
- `-v, --verbose` - Show detailed error messages

**Example:**
```bash
umts show-lineage ./pbmc_mdv/ --packages
```

**Output:**
```
============================================================
LINEAGE INFORMATION: ./pbmc_mdv/
============================================================

UMTS Version: 0.1.0
Created: 2025-11-23T15:30:00Z

Source Files (1):
  1. /path/to/pbmc_seurat.rds
     SHA256: a1b2c3d4e5f6789...
     Size: 45,234,567 bytes
     Modified: 2025-11-20T14:22:00Z

Conversion:
  Function: convert_scanpy_to_mdv
  Parameters:
    max_dims: 3
    delete_existing: True
    label: 
    chunk_data: False
    add_layer_data: True
    gene_identifier_column: None
    n_obs: 2638
    n_vars: 1838

Environment:
  Python: 3.10.12
  Platform: Linux-6.10.14-linuxkit

  Packages:
    mdvtools: 1.0.0
    scanpy: 1.9.3
    scipy: 1.11.2
    numpy: 1.24.3
    pandas: 2.0.3
    anndata: 0.10.0
    h5py: 3.9.0

============================================================
```

## Workflows

### Typical Seurat Workflow

1. **Check dependencies:**
```bash
umts check-deps
```

2. **Convert Seurat file:**
```bash
umts convert-seurat my_seurat.rds ./my_project/
```

3. **Verify lineage:**
```bash
umts show-lineage ./my_project/
```

4. **Serve MDV project:**
```python
from mdvtools import MDVProject
mdv = MDVProject("./my_project/")
mdv.serve()
```

### Batch Processing

Convert multiple files:

```bash
#!/bin/bash
for file in *.rds; do
    name="${file%.rds}"
    echo "Converting $file..."
    umts convert-seurat "$file" "./projects/${name}/"
done
```

### CI/CD Integration

```yaml
# GitHub Actions example
- name: Convert Seurat to MDV
  run: |
    umts check-deps
    umts convert-seurat data/pbmc.rds output/pbmc_mdv/
    
- name: Verify lineage
  run: |
    umts show-lineage output/pbmc_mdv/
```

## Troubleshooting

### "Missing required dependencies" Error

**Problem:** Seurat conversion fails due to missing packages.

**Solution:**
```bash
# Check what's missing
umts check-deps

# Install missing Python packages
pip install rpy2 anndata2ri

# Install R (Ubuntu/Debian)
sudo apt-get install r-base

# Install R packages (in R console)
install.packages("Seurat")
remotes::install_github("mojaveazure/seurat-disk")
```

### "Error converting Seurat file" Error

**Problem:** Conversion fails with error message.

**Solution:**
```bash
# Try with verbose output to see detailed error
umts convert-seurat input.rds output/ --verbose

# Common issues:
# 1. Corrupted RDS file - try loading in R first
# 2. Missing Seurat slots - ensure Seurat object has required data
# 3. Memory issues - try with --chunk-data flag
```

### R Not Found

**Problem:** `check-deps` shows "R not available"

**Solution:**

**Linux:**
```bash
sudo apt-get install r-base
```

**macOS:**
```bash
brew install r
```

**Windows:**
Download from https://www.r-project.org/

### Import Errors

**Problem:** Python cannot find the `umts` command

**Solution:**
```bash
# Reinstall mdvtools
pip install --force-reinstall mdvtools

# Or install in development mode
cd python/
pip install -e .
```

## Help

Get help for any command:

```bash
umts --help
umts convert-seurat --help
umts convert-h5ad --help
umts show-lineage --help
```

## Python API

The CLI tools are also available as Python functions:

```python
from mdvtools.umts import convert_seurat_to_mdv

mdv = convert_seurat_to_mdv(
    seurat_file="input.rds",
    output_folder="./output/",
    max_dims=3,
    track_lineage=True
)
```

See `docs/UMTS_LINEAGE.md` for full Python API documentation.

## Advanced Usage

### Custom Lineage Notes

Add custom notes to lineage after conversion:

```python
from mdvtools import MDVProject
from mdvtools.umts import LineageTracker

# Load existing lineage
mdv = MDVProject("./my_project/")
lineage_data = mdv.get_lineage()

# Add notes
tracker = LineageTracker.from_dict(lineage_data)
tracker.add_note("Applied batch correction")
tracker.add_note("Filtered doublets with Scrublet")
tracker.save("./my_project/")
```

### Verify File Integrity

Check if source file has changed:

```bash
# Get original hash from lineage
umts show-lineage ./my_project/ | grep SHA256

# Compare with current file hash
sha256sum original_file.rds
```

Or in Python:

```python
from mdvtools import MDVProject
from mdvtools.umts import compute_file_hash

mdv = MDVProject("./my_project/")
lineage = mdv.get_lineage()

original_hash = lineage['source_files'][0]['sha256']
current_hash = compute_file_hash('original_file.rds')

if original_hash == current_hash:
    print("✓ File unchanged")
else:
    print("⚠ File has been modified!")
```

## See Also

- [UMTS Lineage Guide](UMTS_LINEAGE.md) - Complete lineage tracking documentation
- [MDVTools Documentation](https://mdv-docs.example.com) - Full MDVTools documentation
- [Tutorial Notebook](../examples/umts_lineage_demo.ipynb) - Interactive tutorial

## Support

For issues or questions:
- GitHub Issues: https://github.com/Taylor-CCB-Group/MDV/issues
- Documentation: docs/UMTS_LINEAGE.md

