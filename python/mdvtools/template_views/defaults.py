"""Hard-coded panels, field-name patterns, and thresholds for template views."""

VIEW_NAMES = [
    "1 Study overview",
    "2 Cell atlas",
    "3 Marker genes and signatures",
    "4 Tissue and disease context",
    "5 RNA-protein concordance",
]

MARKER_PANEL = [
    "CD14",
    "FCGR3A",
    "CD3D",
    "CD3E",
    "MS4A1",
    "CD19",
    "NKG7",
    "GNLY",
    "CD1C",
    "CLEC9A",
    "C1QA",
    "C1QB",
    "LYZ",
    "S100A8",
    "S100A9",
    "VCAN",
    "FCN1",
    "MERTK",
    "LYVE1",
    "APOE",
    "TREM2",
    "IL3RA",
    "KIT",
    "TNF",
    "IL1B",
    "ISG15",
    "CD4",
    "CD8A",
    "FOXP3",
    "MKI67",
]
FEATURE_UMAP_GENES = ["CD14", "FCGR3A", "CD3D", "MS4A1", "CD1C"]
PROTEIN_PANEL = [
    "adt_CD14",
    "CD14",
    "adt_CD16",
    "adt_CD3",
    "adt_CD19",
    "adt_CD4",
    "adt_CD8",
    "adt_CD56",
    "adt_HLA-DR",
    "adt_CD1c",
    "adt_CD11c",
    "adt_CD86",
]

CELL_TYPE_PATTERNS = [
    "annotation_caf",
    "annotation",
    "cell_type",
    "celltype",
    "leiden_res1",
    "leiden",
    "cluster",
    "final_analysis",
]
BROAD_TYPE_PATTERNS = ["sub_bucket_caf", "sub_bucket", "major", "leiden_res1"]
TISSUE_PATTERNS = ["tissue_simple", "tissue", "Run_Tissue_name", "spatial_region"]
DISEASE_PATTERNS = ["diagnosis", "disease"]
TREATMENT_PATTERNS = ["treatment_simple", "treatment"]
RESPONSE_PATTERNS = ["response"]
INFLAMMATION_PATTERNS = ["inflammation_status", "inflammation"]
WORKSTREAM_PATTERNS = ["workstream"]
SEX_PATTERNS = ["gender", "sex"]
SAMPLE_PATTERNS = ["sample_id", "cart_id", "Run_Tissue_name", "slide_ID"]
QC_PATTERNS = [
    "n_genes_by_counts",
    "pct_counts_mt",
    "total_counts",
    "nCount_RNA",
    "nFeature_RNA",
    "nCount_negprobes",
]
FACTOR_PATTERNS = [
    "treatment_simple",
    "treatment",
    "disease_grp_treatment",
    "diagnosis",
    "disease",
    "inflammation_status",
    "inflammation",
    "response",
    "tissue_simple",
    "tissue",
    "Run_Tissue_name",
    "spatial_region",
]
CATEGORICAL_DTYPES = {"text", "text16", "multitext"}
NUMERIC_DTYPES = {"integer", "double", "int32"}
MARKER_DS_NAME = "cluster_markers"
MARKER_CACHE = "cluster_markers.json"
MARKERS_PER_CLUSTER = 20
VARYING_GENES_PER_FACTOR = 15
FEATURE_PLOT_CAP = 24
MISSING_LABELS = {"nd", "nan", "na", "", "undetermined or na", "null"}
SKIP_GENE_PREFIXES = ("MT-", "RPL", "RPS")
SKIP_GENES = {"MALAT1", "NEAT1"}
SEX_GENES = {
    "XIST",
    "DDX3Y",
    "RPS4Y1",
    "EIF1AY",
    "UTY",
    "KDM5D",
    "ZFY",
    "TMSB4Y",
    "NLGN4Y",
    "USP9Y",
}
MATRIX_PREFERENCE = ["rna_logged_counts", "rna_raw_counts", "gs", "rna_expr"]
