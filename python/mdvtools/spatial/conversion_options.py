TABLE_HANDLING_CHOICES = ("merge", "by-region", "per-table")
DEFAULT_TABLE_HANDLING = "by-region"

POINT_TRANSFORM_CHOICES = ("image", "auto", "xenium", "identity", "annotated-element")
DEFAULT_POINT_TRANSFORM = "auto"

# Browser-oriented rewrites of the copied store; see browser_store.py.
BROWSER_PRESET_CHOICES = ("none", "standard", "experimental")
DEFAULT_BROWSER_PRESET = "none"
IMAGE_CODEC_CHOICES = ("none", "htj2k", "jp2k")
IMAGE_PRESET_CHOICES = ("lossless", "balanced", "small")
