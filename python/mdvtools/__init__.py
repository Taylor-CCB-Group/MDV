"""
MDVTools - Multi-Dimensional Viewer Python API

Provides tools for creating, configuring, and managing MDV projects.
"""

from .mdvproject import MDVProject
from .network_helpers import (
    setup_ligand_network,
    create_ligand_network_view,
    setup_cell_communication_network,
)

# UMTS (Universal Modality Translator System) components
try:
    from .umts import LineageTracker, UMTS_VERSION
    _umts_available = True
except ImportError:
    _umts_available = False
    LineageTracker = None
    UMTS_VERSION = None

__all__ = [
    "MDVProject",
    "setup_ligand_network",
    "create_ligand_network_view", 
    "setup_cell_communication_network",
]

# Add UMTS exports if available
if _umts_available:
    __all__.extend([
        "LineageTracker",
        "UMTS_VERSION",
    ])

