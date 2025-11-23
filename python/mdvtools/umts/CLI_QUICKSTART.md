# UMTS CLI Quick Start

Command-line interface for converting single-cell data to MDV projects with automatic lineage tracking.

## Installation

```bash
# Basic (h5ad support only)
pip install mdvtools

# With Seurat support
pip install mdvtools[seurat]
# Also requires R with Seurat and SeuratDisk packages
```

## Quick Examples

### Convert h5ad File

```bash
# Basic conversion
umts convert-h5ad input.h5ad ./output_project/

# With options
umts convert-h5ad input.h5ad ./output_project/ --max-dims 2 --chunk-data
```

### Convert Seurat RDS File

```bash
# Check if Seurat dependencies are installed
umts check-deps

# Convert
umts convert-seurat input.rds ./output_project/

# With options
umts convert-seurat input.rds ./output_project/ --max-dims 3 --label rna_
```

### View Lineage

```bash
# Show provenance information
umts show-lineage ./output_project/

# With package versions
umts show-lineage ./output_project/ --packages
```

### Add External Metadata

```bash
# Add CSV/TSV metadata to existing project
umts add-metadata ./output_project/ phenotypes.csv --key cell_id

# With different join keys
umts add-metadata ./project/ metadata.tsv --key cell_id --metadata-key barcode
```

## Available Commands

```
umts convert-h5ad      # Convert AnnData h5ad file
umts convert-seurat    # Convert Seurat RDS file  
umts add-metadata      # Add CSV/TSV metadata to existing project
umts show-lineage      # Display provenance info
umts check-deps        # Check Seurat dependencies
```

## Common Options

- `--max-dims N` - Maximum dimensions to include (default: 3)
- `--label PREFIX` - Prefix for datasource names
- `--chunk-data` - Lower memory usage (slower)
- `--no-lineage` - Disable lineage tracking
- `--no-delete` - Merge instead of overwrite
- `-v, --verbose` - Show detailed errors

## Help

```bash
umts --help
umts convert-h5ad --help
umts convert-seurat --help
```

## Full Documentation

See [UMTS_CLI.md](../../../../docs/UMTS_CLI.md) for complete documentation.

