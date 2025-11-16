"""
Helper functions for creating network visualizations in MDV.

This module provides convenient functions to set up network charts,
particularly for ligand-receptor interactions and cell-cell communication.
"""

from typing import Optional, List, Dict, Any
import pandas as pd


def setup_ligand_network(
    project,
    interaction_datasource: str,
    cells_datasource: str = "cells",
    ligand_column: str = "ligand_type",
    source_cell_column: str = "source_cell_id",
    target_cell_column: str = "target_cell_id",
    interaction_score_column: str = "interaction_score",
    spatial_distance_column: Optional[str] = None,
    pvalue_column: Optional[str] = None,
    cell_count_column: Optional[str] = None,
    cell_type_column: Optional[str] = None,
    create_view: bool = True,
    view_name: str = "Ligand Network Analysis"
):
    """
    Configure a datasource for ligand-receptor network visualization.
    
    This sets up the required metadata for the CellNetworkChart (Spatial Connectivity Map)
    to display ligand-receptor interactions as a force-directed network.
    
    Parameters
    ----------
    project : MDVProject
        The MDV project instance
    interaction_datasource : str
        Name of the datasource containing interaction data
    cells_datasource : str, default "cells"
        Name of the datasource containing cell/spatial data
    ligand_column : str, default "ligand_type"
        Column containing ligand type/name (used for filtering)
    source_cell_column : str, default "source_cell_id"
        Column containing source cell IDs (cells expressing ligand)
    target_cell_column : str, default "target_cell_id"
        Column containing target cell IDs (cells expressing receptor)
    interaction_score_column : str, default "interaction_score"
        Column with numeric interaction strength (for link thickness/length/color)
    spatial_distance_column : str, optional
        Column with spatial distance between cells (for link length)
    pvalue_column : str, optional
        Column with p-values or FDR (for link color)
    cell_count_column : str, optional
        Column with cell counts (for node size)
    cell_type_column : str, optional
        Column with cell types (for node coloring)
    create_view : bool, default True
        Whether to automatically create a view with the network chart
    view_name : str, default "Ligand Network Analysis"
        Name for the created view
        
    Returns
    -------
    dict
        Configuration dictionary that was set
        
    Examples
    --------
    >>> from mdvtools import MDVProject
    >>> import pandas as pd
    >>> 
    >>> # Create project and load interaction data
    >>> project = MDVProject("my_analysis")
    >>> 
    >>> # Your interaction data
    >>> interactions_df = pd.DataFrame({
    >>>     'ligand_type': ['VEGFA', 'VEGFA', 'TGFB1'],
    >>>     'source_cell_id': ['cell_1', 'cell_1', 'cell_2'],
    >>>     'target_cell_id': ['cell_3', 'cell_4', 'cell_4'],
    >>>     'interaction_score': [0.85, 0.72, 0.91],
    >>>     'spatial_distance': [12.5, 8.3, 15.2],
    >>>     'pvalue': [0.001, 0.003, 0.0001],
    >>>     'cell_count': [15, 15, 22]
    >>> })
    >>> 
    >>> project.add_datasource("ligand_interactions", 
    >>>                         data=interactions_df, 
    >>>                         size=len(interactions_df))
    >>> 
    >>> # Configure for network visualization
    >>> setup_ligand_network(
    >>>     project,
    >>>     interaction_datasource="ligand_interactions",
    >>>     ligand_column="ligand_type",
    >>>     source_cell_column="source_cell_id",
    >>>     target_cell_column="target_cell_id",
    >>>     interaction_score_column="interaction_score",
    >>>     spatial_distance_column="spatial_distance",
    >>>     pvalue_column="pvalue",
    >>>     cell_count_column="cell_count"
    >>> )
    >>> 
    >>> # Now in the GUI: Add Chart → Select "ligand_interactions" datasource
    >>> #                → "Spatial Connectivity Map" will appear
    >>> #                → Select ligand type → Network displays!
    """
    
    # Determine which columns to use for visualization
    link_length_col = spatial_distance_column or interaction_score_column
    link_thickness_col = interaction_score_column
    link_color_col = pvalue_column or interaction_score_column
    node_size_col = cell_count_column or interaction_score_column
    
    # Check all required columns exist
    required_cols = [
        ligand_column,
        source_cell_column,
        target_cell_column,
        link_length_col,
        link_thickness_col,
        link_color_col,
        node_size_col
    ]
    
    if cell_type_column:
        required_cols.append(cell_type_column)
    
    missing_cols = project.check_columns_exist(interaction_datasource, required_cols)
    if missing_cols:
        raise ValueError(
            f"Missing columns in '{interaction_datasource}': {', '.join(missing_cols)}"
        )
    
    # Configure the interactions metadata
    md = project.get_datasource_metadata(interaction_datasource)
    md["interactions"] = {
        "pivot_column": ligand_column,
        "is_single_region": False,  # Network can span multiple regions
        "interaction_columns": [source_cell_column, target_cell_column],
        "spatial_connectivity_map": {
            "link_length": link_length_col,
            "link_thickness": link_thickness_col,
            "link_color": link_color_col,
            "node_size": node_size_col,
        },
    }
    
    # Add cell type if specified
    if cell_type_column:
        md["interactions"]["spatial_connectivity_map"]["cell_type"] = cell_type_column
    
    project.set_datasource_metadata(md)
    
    # Set up linking to cells datasource if it exists
    try:
        parent_md = project.get_datasource_metadata(cells_datasource)
        project.insert_link(
            interaction_datasource,
            cells_datasource,
            "interactions",
            {
                "interaction_columns": [source_cell_column, target_cell_column],
                "pivot_column": ligand_column,
                "is_single_region": False,
            },
        )
    except AttributeError:
        print(f"Warning: '{cells_datasource}' datasource not found. Skipping link creation.")
    
    # Create view with network chart
    if create_view:
        create_ligand_network_view(
            project,
            interaction_datasource,
            cells_datasource,
            view_name
        )
    
    print(f"✓ Network configuration complete for '{interaction_datasource}'")
    print(f"  Pivot: {ligand_column}")
    print(f"  Nodes: {source_cell_column} ↔ {target_cell_column}")
    print(f"  Visual encoding:")
    print(f"    - Link length: {link_length_col}")
    print(f"    - Link thickness: {link_thickness_col}")
    print(f"    - Link color: {link_color_col}")
    print(f"    - Node size: {node_size_col}")
    if cell_type_column:
        print(f"    - Node color: {cell_type_column}")
    
    return md["interactions"]


def create_ligand_network_view(
    project,
    interaction_datasource: str,
    cells_datasource: str = "cells",
    view_name: str = "Ligand Network Analysis",
    include_spatial_panel: bool = True,
    network_panel_width: int = 60,
    spatial_panel_width: int = 40
):
    """
    Create a view with network chart and optional spatial centroid plots.
    
    Creates a split-view layout with:
    - Left panel: Network visualization of interactions
    - Right panel: Spatial centroid plots (optional)
    
    Parameters
    ----------
    project : MDVProject
        The MDV project instance
    interaction_datasource : str
        Name of the datasource with interaction metadata configured
    cells_datasource : str, default "cells"
        Name of the datasource containing spatial cell data
    view_name : str, default "Ligand Network Analysis"
        Name for the view
    include_spatial_panel : bool, default True
        Whether to include spatial visualization panel
    network_panel_width : int, default 60
        Width percentage for network panel (0-100)
    spatial_panel_width : int, default 40
        Width percentage for spatial panel (0-100)
        
    Returns
    -------
    None
    
    Examples
    --------
    >>> # After setting up network with setup_ligand_network()
    >>> create_ligand_network_view(
    >>>     project,
    >>>     "ligand_interactions",
    >>>     "cells",
    >>>     view_name="VEGFA Signaling Network"
    >>> )
    """
    
    # Get interaction metadata to find ligand types
    md = project.get_datasource_metadata(interaction_datasource)
    interactions = md.get("interactions")
    
    if not interactions:
        raise ValueError(
            f"Datasource '{interaction_datasource}' does not have interaction metadata. "
            "Run setup_ligand_network() first."
        )
    
    # Get ligand types for selection dialog
    pivot_column = interactions["pivot_column"]
    ligand_col_md = project.get_column_metadata(interaction_datasource, pivot_column)
    ligand_types = ligand_col_md.get("values", [])
    
    # Create network panel charts
    network_charts = []
    
    # Add selection dialog for filtering
    filter_dialog = project.get_selection_dialog(
        interaction_datasource,
        [
            {"column": pivot_column},
            {"column": interactions["interaction_columns"][0]},
            {"column": interactions["interaction_columns"][1]},
        ]
    )
    filter_dialog["position"] = [5, 5]
    filter_dialog["size"] = [300, 250]
    network_charts.append(filter_dialog)
    
    # Note: CellNetworkChart must be added manually via GUI
    # because it requires selecting a specific ligand type
    # We provide instructions in the view
    info_text = {
        "type": "text_box_chart",
        "param": [],
        "text": (
            f"<h2>Ligand-Receptor Network</h2>\n"
            f"<p><b>To add network chart:</b></p>\n"
            f"<ol>\n"
            f"  <li>Click 'Add Chart' button</li>\n"
            f"  <li>Select datasource: <code>{interaction_datasource}</code></li>\n"
            f"  <li>Chart type: <b>Spatial Connectivity Map</b></li>\n"
            f"  <li>Select ligand type: {', '.join(ligand_types[:5])}"
            f"{'...' if len(ligand_types) > 5 else ''}</li>\n"
            f"  <li>Click 'Add Chart'</li>\n"
            f"</ol>\n"
            f"<p><b>Configured visual encoding:</b></p>\n"
            f"<ul>\n"
            f"  <li>Link length: {interactions['spatial_connectivity_map']['link_length']}</li>\n"
            f"  <li>Link thickness: {interactions['spatial_connectivity_map']['link_thickness']}</li>\n"
            f"  <li>Link color: {interactions['spatial_connectivity_map']['link_color']}</li>\n"
            f"  <li>Node size: {interactions['spatial_connectivity_map']['node_size']}</li>\n"
            f"</ul>"
        ),
        "position": [310, 5],
        "size": [500, 350]
    }
    network_charts.append(info_text)
    
    # Create spatial panel (if cells datasource exists and option enabled)
    spatial_charts = []
    if include_spatial_panel:
        try:
            cells_md = project.get_datasource_metadata(cells_datasource)
            # Add placeholder text for spatial visualizations
            spatial_info = {
                "type": "text_box_chart",
                "param": [],
                "text": (
                    f"<h2>Spatial Context</h2>\n"
                    f"<p>Add centroid plots here to see spatial distribution of interacting cells.</p>\n"
                    f"<p>Click 'Add Chart' → Select '{cells_datasource}' → Choose spatial region.</p>"
                ),
                "position": [5, 5],
                "size": [400, 200]
            }
            spatial_charts.append(spatial_info)
        except AttributeError:
            print(f"Warning: No '{cells_datasource}' datasource found. Skipping spatial panel.")
            include_spatial_panel = False
    
    # Configure view
    view_config = {
        "initialCharts": {
            interaction_datasource: network_charts,
        },
        "dataSources": {
            interaction_datasource: {"panelWidth": 100 if not include_spatial_panel else network_panel_width},
        }
    }
    
    if include_spatial_panel and spatial_charts:
        view_config["initialCharts"][cells_datasource] = spatial_charts
        view_config["dataSources"][cells_datasource] = {"panelWidth": spatial_panel_width}
    
    project.set_view(view_name, view_config)
    
    print(f"✓ View '{view_name}' created")
    print(f"  Network panel: {network_panel_width}% width")
    if include_spatial_panel:
        print(f"  Spatial panel: {spatial_panel_width}% width")


def setup_cell_communication_network(
    project,
    interaction_datasource: str,
    communication_type_column: str = "communication_type",
    sender_column: str = "sender_cell",
    receiver_column: str = "receiver_cell",
    score_column: str = "communication_score",
    **kwargs
):
    """
    Convenience wrapper for setup_ligand_network with cell communication naming.
    
    This is identical to setup_ligand_network but uses more general terminology
    for cell-cell communication studies.
    
    Parameters
    ----------
    project : MDVProject
        The MDV project instance
    interaction_datasource : str
        Name of the datasource containing communication data
    communication_type_column : str, default "communication_type"
        Column containing communication pathway/type
    sender_column : str, default "sender_cell"
        Column containing sender cell IDs
    receiver_column : str, default "receiver_cell"
        Column containing receiver cell IDs
    score_column : str, default "communication_score"
        Column with communication strength score
    **kwargs
        Additional arguments passed to setup_ligand_network()
        
    Examples
    --------
    >>> setup_cell_communication_network(
    >>>     project,
    >>>     "cell_comm",
    >>>     communication_type_column="pathway",
    >>>     sender_column="from_cell",
    >>>     receiver_column="to_cell",
    >>>     score_column="score"
    >>> )
    """
    return setup_ligand_network(
        project,
        interaction_datasource,
        ligand_column=communication_type_column,
        source_cell_column=sender_column,
        target_cell_column=receiver_column,
        interaction_score_column=score_column,
        **kwargs
    )

