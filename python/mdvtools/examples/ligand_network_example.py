"""
Example: Creating a Ligand-Receptor Network Visualization in MDV

This example demonstrates how to:
1. Prepare ligand-receptor interaction data
2. Configure it for network visualization
3. Create an interactive view with the network chart
4. Link to spatial data
"""

import pandas as pd
import numpy as np
from mdvtools import MDVProject
from mdvtools.network_helpers import setup_ligand_network

# =============================================================================
# Step 1: Prepare Example Data
# =============================================================================

# Create sample ligand-receptor interaction data
np.random.seed(42)

# Define some ligands, receptors, and cells
ligands = ['VEGFA', 'VEGFA', 'VEGFA', 'TGFB1', 'TGFB1', 'TGFB1', 'PDGFB', 'PDGFB', 'FGF2', 'FGF2']
receptors = ['KDR', 'FLT1', 'FLT1', 'TGFBR1', 'TGFBR2', 'TGFBR1', 'PDGFRB', 'PDGFRA', 'FGFR1', 'FGFR2']
source_cells = ['T_cell_1', 'T_cell_1', 'T_cell_2', 'Cancer_A', 'Cancer_A', 'Cancer_B', 
                'Fibroblast_1', 'Fibroblast_1', 'Endothelial_1', 'Endothelial_1']
target_cells = ['Endothelial_1', 'Endothelial_2', 'Endothelial_1', 'Fibroblast_1', 'Fibroblast_2',
                'Fibroblast_1', 'Cancer_A', 'Cancer_B', 'T_cell_1', 'T_cell_2']
cell_types_source = ['T_cell', 'T_cell', 'T_cell', 'Cancer', 'Cancer', 'Cancer', 
                      'Fibroblast', 'Fibroblast', 'Endothelial', 'Endothelial']
cell_types_target = ['Endothelial', 'Endothelial', 'Endothelial', 'Fibroblast', 'Fibroblast',
                      'Fibroblast', 'Cancer', 'Cancer', 'T_cell', 'T_cell']

# Create interaction dataframe
interactions_df = pd.DataFrame({
    'ligand_type': ligands,
    'receptor_type': receptors,
    'source_cell_id': source_cells,
    'target_cell_id': target_cells,
    'source_cell_type': cell_types_source,
    'target_cell_type': cell_types_target,
    'interaction_score': np.random.uniform(0.5, 1.0, 10),
    'spatial_distance': np.random.uniform(5, 50, 10),
    'pvalue': np.random.uniform(0.0001, 0.05, 10),
    'fdr': np.random.uniform(0.001, 0.1, 10),
    'cell_count': np.random.randint(5, 30, 10),
})

print("Sample interaction data:")
print(interactions_df.head())
print(f"\nTotal interactions: {len(interactions_df)}")
print(f"Unique ligands: {interactions_df['ligand_type'].nunique()}")
print(f"Unique cells: {len(set(source_cells) | set(target_cells))}")

# =============================================================================
# Step 2: Create MDV Project and Add Data
# =============================================================================

project = MDVProject(
    "ligand_network_example",
    delete_existing=True  # Start fresh
)

# Add interaction datasource
project.add_datasource(
    "ligand_interactions",
    dataframe=interactions_df
)

print(f"\n✓ Added datasource 'ligand_interactions' with {len(interactions_df)} rows")

# Optional: Add cell spatial data if you have it
# This enables linking between network and spatial views
has_spatial_data = False

try:
    # Example spatial data (normally you'd load from actual data)
    all_cells = list(set(source_cells) | set(target_cells))
    cells_df = pd.DataFrame({
        'cell_id': all_cells,
        'x': np.random.uniform(0, 1000, len(all_cells)),
        'y': np.random.uniform(0, 1000, len(all_cells)),
        'cell_type': [
            'T_cell' if 'T_cell' in c else
            'Cancer' if 'Cancer' in c else
            'Fibroblast' if 'Fibroblast' in c else
            'Endothelial'
            for c in all_cells
        ]
    })
    
    project.add_datasource("cells", dataframe=cells_df)
    
    # Set up spatial regions (required for spatial charts)
    project.set_region_data(
        "cells",
        {"roi_1": {"width": 1000, "height": 1000}},
        region_field="cell_id",  # All cells in same region for this example
        default_color="cell_type",
        position_fields=["x", "y"],
        scale_unit="µm",
        scale=1.0
    )
    
    has_spatial_data = True
    print("✓ Added spatial data for cells")
    
except Exception as e:
    print(f"Note: Skipping spatial data ({e})")

# =============================================================================
# Step 3: Configure for Network Visualization
# =============================================================================

print("\n" + "="*60)
print("Configuring network visualization...")
print("="*60)

config = setup_ligand_network(
    project,
    interaction_datasource="ligand_interactions",
    cells_datasource="cells" if has_spatial_data else None,
    ligand_column="ligand_type",
    source_cell_column="source_cell_id",
    target_cell_column="target_cell_id",
    interaction_score_column="interaction_score",
    spatial_distance_column="spatial_distance",
    pvalue_column="pvalue",
    cell_count_column="cell_count",
    cell_type_column="source_cell_type",
    create_view=True,
    view_name="Ligand-Receptor Network"
)

# =============================================================================
# Step 4: Save Project
# =============================================================================

project.save()

print("\n" + "="*60)
print("Project created successfully!")
print("="*60)
print(f"\nProject location: {project.dir}")
print("\nNext steps:")
print("1. Open MDV GUI and load this project")
print("2. Navigate to 'Ligand-Receptor Network' view")
print("3. Click 'Add Chart' button")
print("4. Select datasource: 'ligand_interactions'")
print("5. Choose chart type: 'Spatial Connectivity Map'")
print("6. Select a ligand type (e.g., 'VEGFA', 'TGFB1', or 'PDGFB')")
print("7. Click 'Add Chart' - the network will appear!")
print("\n" + "="*60)
print("Chart Features:")
print("="*60)
print("- Nodes represent cells (source/target)")
print("- Links represent ligand-receptor interactions")
print("- Link thickness shows interaction strength")
print("- Link color shows statistical significance (p-value)")
print("- Node size shows number of interactions")
print("- Drag nodes to rearrange")
print("- Click links to highlight in other charts")
print("- Use settings panel to adjust:")
print("  * Link strength (physics)")
print("  * Node repulsion (spacing)")
print("  * Visual scales (thickness, color, size)")
print("  * Center cell (ego network mode)")

# =============================================================================
# Optional: Advanced Usage
# =============================================================================

def create_multiple_network_views(project):
    """
    Create separate views for each ligand type.
    """
    ligand_types = interactions_df['ligand_type'].unique()
    
    for ligand in ligand_types:
        view_config = {
            "initialCharts": {
                "ligand_interactions": [
                    {
                        "type": "text_box_chart",
                        "param": [],
                        "text": f"<h2>{ligand} Network</h2><p>Add 'Spatial Connectivity Map' chart and select '{ligand}'</p>",
                        "position": [5, 5],
                        "size": [400, 200]
                    }
                ]
            }
        }
        
        project.set_view(f"{ligand} Network", view_config)
        print(f"✓ Created view for {ligand}")

# Uncomment to create per-ligand views:
# create_multiple_network_views(project)
# project.save()

# =============================================================================
# Alternative: Using CellTalkDB or similar database results
# =============================================================================

def load_celltalkdb_results(celltalkdb_file: str):
    """
    Example of loading pre-computed cell-cell communication results.
    
    CellTalkDB, CellPhoneDB, CellChat, NicheNet, etc. typically output:
    - ligand/receptor pairs
    - cell type pairs
    - interaction scores
    - p-values
    
    Adapt column names to match your tool's output format.
    """
    # Example structure
    df = pd.read_csv(celltalkdb_file)
    
    # Map to MDV expected columns
    df_mapped = df.rename(columns={
        'ligand': 'ligand_type',
        'receptor': 'receptor_type',
        'cell_type_a': 'source_cell_id',
        'cell_type_b': 'target_cell_id',
        'score': 'interaction_score',
        'p_value': 'pvalue'
    })
    
    return df_mapped

# Example usage:
# df = load_celltalkdb_results("celltalkdb_output.csv")
# project.add_datasource("cell_communication", data=df, size=len(df))
# setup_ligand_network(project, "cell_communication", ...)

