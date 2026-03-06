# MDV Differential Gene Expression (DGE) Implementation

This document describes the "Tier 1" browser-based Differential Gene Expression (DGE) implementation in MDV.

## 1. Overview

The MDV DGE engine (`mdv-dge`) provides high-performance, real-time statistical analysis of single-cell RNA-seq data directly in the browser. It is designed to work with hundreds of thousands of cells and tens of thousands of genes without requiring a dedicated high-memory backend for computation.

## 2. Architecture

The DGE feature is split into four layers:

1.  **UI Layer (`DGEDialogReact.tsx`)**: A React/MUI dialog that allows users to select a grouping column (e.g., Leiden cluster, Disease state) and compare a target group against a reference group (or "rest").
2.  **Integration Layer (`dgeIntegration.ts`)**: The "glue" that connects MDV's `DataStore` and `ChartManager` to the DGE engine. It handles data discovery via `rows_as_columns` links, coordinates batch loading of gene data, and writes results back as new columns.
3.  **Calculation Engine (`DGEDimension.ts` & `dgeWorker.ts`)**: Orchestrates the DGE pipeline. It divides genes into batches and can optionally offload calculations to **Web Workers** to keep the UI responsive.
4.  **Statistical Core (`dgeStats.ts`)**: Pure TypeScript implementation of the statistical algorithms (Welch's t-test, Welford's online variance, Benjamini-Hochberg correction).

## 3. Statistical Methodology

### Welch's t-test
We use **Welch’s t-test** (unequal variances t-test) to compare gene expression between two groups. This is the same default method used by popular tools like **Scanpy**.

*   **P-value Calculation**: Computed using the Student's t-distribution CDF (via the Regularized Incomplete Beta function).
*   **Multiple Testing Correction**: P-values are adjusted using the **Benjamini-Hochberg (BH)** procedure to control the False Discovery Rate (FDR).

### Online Variance (Welford's Algorithm)
To handle large datasets without multiple passes through the data, we use **Welford’s online algorithm**. This allows us to calculate the mean and variance of gene expression in a single pass over the cell indices, minimizing memory pressure.

### Adaptive Effect Size
The engine automatically detects the normalization of the input data and selects the appropriate effect size metric:

*   **Log1p-normalized data**: Uses **Log2 Fold Change (Log2FC)** with `expm1` back-transformation.
*   **Linear / raw count data**: Uses **Log2 Fold Change (Log2FC)** directly on means (no `expm1`).
*   **Z-scored data**: Uses **Mean Difference** (`meanTarget - meanReference`).

Detection is automatic via a multi-gene probe (see gotcha below).

### Gotcha: Data Type Detection for Sparse Raw Counts

> **If the volcano plot shows many genes clamped at exactly ±10 effect size, or highly-expressed housekeeping genes have implausible log2FC values, check that data type detection is working correctly.**

The engine probes the first 50 genes to decide if expression values are log1p-normalized (small non-negative, max ≤ 20), linear/raw counts (non-negative, max > 20), or z-scored (has negatives). For **sparse raw count data**, most genes have very low counts (1–3 per cell) which look identical to log1p-normalized values. If detection breaks at the first gene it finds, it may misclassify raw counts as log1p — causing `expm1()` to be applied to already-linear values, which exponentially inflates fold changes for highly-expressed genes.

**The fix** (implemented in `dgeIntegration.ts`): the probe loop scans across **all** 50 sampled genes, tracking the global maximum value, and only decides after examining all of them. It exits early only when it finds a definitive signal (a negative value → z-scored, or a value > 20 → linear).

Even with correct detection, genes expressed exclusively in one group (e.g., 72 UC cells, 0 Healthy cells) will produce mathematically extreme fold changes that are clamped to ±10. This is expected for very sparse genes and can be mitigated with a minimum cell count filter if needed.

## 4. Performance & Scalability

### Batching
Data is requested from the backend in batches (default: **2000 genes per batch**). This strikes a balance between network overhead and browser memory usage.

### Memory Optimization
*   **SharedArrayBuffer**: Expression data is loaded into `SharedArrayBuffer` objects, allowing efficient sharing between the main thread and Web Workers.
*   **Sparse Data Handling**: In sparse datasets, `NaN` values in MDV's dense buffers are treated as `0` during statistical accumulation, ensuring correct cell counts for variance calculations.

### Browser Limits
*   **Memory**: Each gene for 60k cells takes ~240KB. A batch of 2000 genes takes ~480MB.
*   **Dataset Size**: While Tier 1 works well for up to 500k cells, datasets of **10M+ cells** require either sub-sampling or server-side execution (Tier 2) due to browser memory and transfer limits.

## 5. Backend Integration: HTTP Streaming

To prevent the backend from crashing when requesting thousands of columns at once, we implemented **HTTP Streaming** (Chunked Transfer Encoding) in the Flask server:

*   **`yield_byte_data`**: The Python backend streams individual gene columns from HDF5 as they are read, instead of buffering the entire multi-gigabyte response in RAM.
*   **No Content-Length**: This removes the risk of `ERR_CONTENT_LENGTH_MISMATCH` errors caused by large binary payloads.

## 6. Configuration & Troubleshooting

### Adjusting Batch Size
If the UI feels sluggish or you encounter memory errors, you can adjust the default DGE batch size in the following centralized location:

*   **File:** `src/lib/constants.ts`
*   **Constant:** `DEFAULT_DGE_BATCH_SIZE`

This single source of truth is used as the default fallback in both the integration layer (`dgeIntegration.ts`) and the core DGE engine (`DGEDimension.ts`).

**After making changes** you must perform a hard browser refresh (Cmd+Shift+R or Ctrl+F5) to clear the cached JavaScript.

#### Guidance on choosing a value

| Batch Size | Effect |
| :--- | :--- |
| **500** | Safer for large datasets (>100k cells). Lower peak memory per batch (~120MB). More HTTP round-trips. |
| **1000** | Good balance for datasets up to ~200k cells. |
| **2000** *(default)* | Optimal for datasets up to ~60k cells. Reduces round-trips significantly. |
| **5000+** | Only suitable for small datasets (<20k cells). Risk of browser memory pressure / GC pauses. |

### Visualizing Results
After a DGE run completes, a **Volcano Plot** is automatically generated. You can also manually create scatter plots or color other charts using the newly generated `dge_effect_size` and `dge_neg_log10_pval_adj` columns.

### Caching
DGE results are cached in memory. If you run DGE a second time on the same dataset, the browser will reuse the cached gene expression data, making the second run nearly instantaneous.
