TABLE_HANDLING_CHOICES = ("merge", "by-region", "per-table")
DEFAULT_TABLE_HANDLING = "by-region"

POINT_TRANSFORM_CHOICES = ("image", "auto", "xenium", "identity", "annotated-element")
DEFAULT_POINT_TRANSFORM = "auto"

# Viv's tile size is the chunk size, and Xenium stores written by spatialdata_io use 4096 x 4096.
DEFAULT_IMAGE_CHUNK_SIZE = 1024
