"""
UMTS (Universal Modality Translator System) for MDVTools.

Provides lineage tracking and provenance recording for data conversions.

Phase 1 Implementation: Lineage Tracking
- Track source files with SHA256 hashes
- Record conversion parameters and environment
- Maintain audit trail for reproducibility

Future phases will add:
- Global ID system for entity tracking
- Schema governance framework
- Translator registry for extensibility
- AI/MCP adaptive merging
"""

from .lineage import LineageTracker, UMTS_VERSION
from .utils import (
    compute_file_hash,
    get_package_versions,
    get_timestamp,
    get_file_metadata,
    get_environment_info
)

# Optional Seurat support (requires rpy2 and anndata2ri)
try:
    from .seurat_reader import convert_seurat_to_mdv, check_seurat_dependencies
    _seurat_available = True
except ImportError:
    _seurat_available = False
    convert_seurat_to_mdv = None
    check_seurat_dependencies = None

__version__ = UMTS_VERSION

__all__ = [
    'LineageTracker',
    'UMTS_VERSION',
    'compute_file_hash',
    'get_package_versions',
    'get_timestamp',
    'get_file_metadata',
    'get_environment_info',
]

# Add Seurat exports if available
if _seurat_available:
    __all__.extend([
        'convert_seurat_to_mdv',
        'check_seurat_dependencies',
    ])

