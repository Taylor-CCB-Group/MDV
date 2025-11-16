"""
Flexible Network Chart - Zero-Friction Example

This demonstrates how easy it is to visualize networks WITHOUT any Python configuration!
Just load your data and select columns in the GUI.
"""

import pandas as pd
import numpy as np
from mdvtools import MDVProject

# =============================================================================
# Create Sample Network Data (Any CSV format works!)
# =============================================================================

np.random.seed(42)

# Example 1: Simple interaction network
simple_network = pd.DataFrame({
    'from_node': ['A', 'A', 'B', 'B', 'C', 'D', 'D', 'E'],
    'to_node': ['B', 'C', 'C', 'D', 'D', 'E', 'F', 'F'],
    'weight': np.random.uniform(0.5, 1.0, 8),
})

# Example 2: Ligand-receptor interactions (no special format!)
ligand_receptor = pd.DataFrame({
    'ligand_cell': ['T_cell_1', 'T_cell_1', 'T_cell_2', 'Cancer_A', 'Cancer_A', 
                    'Fibroblast_1', 'Fibroblast_1', 'Endothelial_1'],
    'receptor_cell': ['Endothelial_1', 'Endothelial_2', 'Endothelial_1', 
                      'Fibroblast_1', 'Fibroblast_2', 'Cancer_A', 'Cancer_B', 'T_cell_1'],
    'interaction_score': np.random.uniform(0.6, 1.0, 8),
    'ligand_type': ['VEGFA', 'VEGFA', 'VEGFA', 'TGFB1', 'TGFB1', 'PDGFB', 'PDGFB', 'FGF2'],
    'spatial_distance': np.random.uniform(5, 50, 8),
    'pvalue': np.random.uniform(0.0001, 0.05, 8),
    'cell_count': np.random.randint(5, 30, 8),
    'cell_type': ['T_cell', 'T_cell', 'T_cell', 'Cancer', 'Cancer', 
                  'Fibroblast', 'Fibroblast', 'Endothelial'],
})

# Example 3: Gene regulatory network
gene_network = pd.DataFrame({
    'transcription_factor': ['TP53', 'TP53', 'STAT3', 'STAT3', 'NF-KB', 'NF-KB'],
    'target_gene': ['MDM2', 'BAX', 'IL6', 'SOCS3', 'IL1B', 'TNF'],
    'regulatory_score': [0.95, 0.88, 0.92, 0.87, 0.90, 0.85],
    'regulation_type': ['activation', 'activation', 'activation', 'activation', 'activation', 'activation'],
})

print("="*60)
print("FLEXIBLE NETWORK CHART - Example Data")
print("="*60)
print("\n1. Simple Network:")
print(simple_network)
print(f"\nNodes: {len(set(simple_network['from_node']) | set(simple_network['to_node']))}")
print(f"Edges: {len(simple_network)}")

print("\n2. Ligand-Receptor Network:")
print(ligand_receptor.head())
print(f"\nNodes: {len(set(ligand_receptor['ligand_cell']) | set(ligand_receptor['receptor_cell']))}")
print(f"Edges: {len(ligand_receptor)}")
print(f"Ligand types: {ligand_receptor['ligand_type'].unique()}")

print("\n3. Gene Regulatory Network:")
print(gene_network)

# =============================================================================
# Load into MDV - Just add datasources, NO CONFIGURATION!
# =============================================================================

project = MDVProject("flexible_network_examples", delete_existing=True)

# Example 1: Simple network
project.add_datasource("simple", dataframe=simple_network)

# Example 2: Ligand-receptor (all columns included)
project.add_datasource("ligand_receptor", dataframe=ligand_receptor)

# Example 3: Gene network
project.add_datasource("genes", dataframe=gene_network)

# Create views with instructions
project.set_view("Simple Network", {
    "initialCharts": {
        "simple": [{
            "type": "text_box_chart",
            "param": [],
            "text": """
<h2>Simple Network Example</h2>
<p><b>To visualize:</b></p>
<ol>
  <li>Click <b>Add Chart</b></li>
  <li>Select datasource: <code>simple</code></li>
  <li>Chart type: <b>Network Graph (Flexible)</b></li>
  <li>Select columns:
    <ul>
      <li>Source Node: <code>from_node</code></li>
      <li>Target Node: <code>to_node</code></li>
      <li>Edge Weight: <code>weight</code></li>
      <li>Leave optional columns empty for now</li>
    </ul>
  </li>
  <li>Click <b>Add Chart</b> - Done!</li>
</ol>
<p><b>That's it!</b> No Python configuration needed.</p>
<p>Drag nodes to rearrange. Click edges to highlight.</p>
            """,
            "position": [5, 5],
            "size": [500, 400]
        }]
    }
})

project.set_view("Ligand-Receptor Network", {
    "initialCharts": {
        "ligand_receptor": [{
            "type": "text_box_chart",
            "param": [],
            "text": """
<h2>Ligand-Receptor Network Example</h2>
<p><b>Full configuration with ALL optional columns:</b></p>
<ol>
  <li>Click <b>Add Chart</b></li>
  <li>Select datasource: <code>ligand_receptor</code></li>
  <li>Chart type: <b>Network Graph (Flexible)</b></li>
  <li>Select columns:
    <ul>
      <li><b>Source Node:</b> <code>ligand_cell</code></li>
      <li><b>Target Node:</b> <code>receptor_cell</code></li>
      <li><b>Edge Weight:</b> <code>interaction_score</code></li>
      <li><b>Category Filter:</b> <code>ligand_type</code></li>
      <li><b>Edge Length:</b> <code>spatial_distance</code></li>
      <li><b>Edge Color:</b> <code>pvalue</code></li>
      <li><b>Node Size:</b> <code>cell_count</code></li>
      <li><b>Node Type:</b> <code>cell_type</code></li>
    </ul>
  </li>
  <li>Click <b>Add Chart</b></li>
</ol>
<p><b>Result:</b> Fully configured network with:</p>
<ul>
  <li>Link thickness = interaction strength</li>
  <li>Link length = spatial distance</li>
  <li>Link color = p-value (red=significant)</li>
  <li>Node size = cell count</li>
  <li>Node color = cell type</li>
</ul>
            """,
            "position": [5, 5],
            "size": [550, 500]
        }]
    }
})

project.set_view("Gene Regulatory Network", {
    "initialCharts": {
        "genes": [{
            "type": "text_box_chart",
            "param": [],
            "text": """
<h2>Gene Regulatory Network Example</h2>
<p><b>Minimal configuration (3 required columns only):</b></p>
<ol>
  <li>Click <b>Add Chart</b></li>
  <li>Select datasource: <code>genes</code></li>
  <li>Chart type: <b>Network Graph (Flexible)</b></li>
  <li>Select columns:
    <ul>
      <li>Source Node: <code>transcription_factor</code></li>
      <li>Target Node: <code>target_gene</code></li>
      <li>Edge Weight: <code>regulatory_score</code></li>
      <li>(Optional) Category: <code>regulation_type</code></li>
    </ul>
  </li>
  <li>Click <b>Add Chart</b></li>
</ol>
<p><b>Tip:</b> Enable "Show Directionality" in settings to see TF → gene arrows.</p>
            """,
            "position": [5, 5],
            "size": [500, 400]
        }]
    }
})

project.save()

print("\n" + "="*60)
print("PROJECT CREATED!")
print("="*60)
print(f"\nLocation: {project.dir}")
print("\nDatasets loaded:")
print("  1. 'simple' - Basic network (3 columns)")
print("  2. 'ligand_receptor' - Full network (8 columns)")  
print("  3. 'genes' - Gene regulatory network (4 columns)")

print("\n" + "="*60)
print("NEXT STEPS - Open in MDV GUI:")
print("="*60)
print("\n1. Launch MDV and load this project")
print("2. Navigate to any view (Simple/Ligand-Receptor/Gene)")
print("3. Click 'Add Chart' button")
print("4. Select datasource from dropdown")
print("5. Choose chart type: 'Network Graph (Flexible)'")
print("6. Select columns via dropdowns:")
print("   - Required: Source Node, Target Node, Edge Weight")
print("   - Optional: 5 more columns for advanced features")
print("7. Click 'Add Chart' - Network appears instantly!")

print("\n" + "="*60)
print("KEY DIFFERENCES FROM METADATA-DRIVEN CHART:")
print("="*60)
print("✅ NO Python configuration required")
print("✅ Works with ANY table format")
print("✅ Select columns directly in GUI")
print("✅ No metadata, no setup, no friction")
print("✅ Perfect for exploration and varied data")

print("\n" + "="*60)
print("COMPARISON:")
print("="*60)
print("\nMetadata-driven (CellNetworkChart):")
print("  - Requires: Python setup with set_interactions()")
print("  - Datasource: Must have 'interactions' metadata")
print("  - GUI: 1 dropdown (select category)")
print("  - Use case: Production, standardized pipelines")

print("\nFlexible (This chart):")
print("  - Requires: Nothing! Just data")
print("  - Datasource: Any table with network data")
print("  - GUI: 8 dropdowns (select all columns)")
print("  - Use case: Exploration, varied formats")

print("\n" + "="*60)
print("READY TO VISUALIZE!")
print("="*60)
print("\nBoth charts are available - choose based on your workflow:")
print("  - Quick exploration? → Use Flexible Network")
print("  - Production pipeline? → Use Spatial Connectivity Map")

