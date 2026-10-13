"""
MDVTools - Multi-Dimensional Viewer Python API

Provides tools for creating, configuring, and managing MDV projects.
"""

from importlib import import_module
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .mdvproject import MDVProject
    from .network_helpers import (
        setup_ligand_network,
        create_ligand_network_view,
        setup_cell_communication_network,
    )

# Exported names are imported on first access, so importing a submodule such as
# mdvtools.jobs.run_worker does not load MDVProject and its stack (pandas, scanpy).
_LAZY_EXPORTS = {
    "MDVProject": ".mdvproject",
    "setup_ligand_network": ".network_helpers",
    "create_ligand_network_view": ".network_helpers",
    "setup_cell_communication_network": ".network_helpers",
}

__all__ = list(_LAZY_EXPORTS)


def __getattr__(name: str):
    if name not in _LAZY_EXPORTS:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    value = getattr(import_module(_LAZY_EXPORTS[name], __name__), name)
    globals()[name] = value
    return value
