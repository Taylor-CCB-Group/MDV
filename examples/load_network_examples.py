"""
Quick script to load example network CSV files into MDV

Run this to create a project with all 3 example networks ready to visualize!
"""

import pandas as pd
import os
from pathlib import Path

# Add mdvtools to path if running from examples directory
import sys
sys.path.insert(0, str(Path(__file__).parent.parent / "python"))

from mdvtools import MDVProject

# Get the data directory
data_dir = Path(__file__).parent / "data"

print("="*70)
print("LOADING NETWORK EXAMPLE DATA")
print("="*70)

# Create project
project = MDVProject(
    str(Path(__file__).parent / "network_examples"),
    delete_existing=True
)

# Load all 3 CSV files
examples = [
    {
        "file": "ligand_receptor_network_example.csv",
        "name": "ligand_receptor",
        "description": "Full-featured: 30 interactions with all visual encodings"
    },
    {
        "file": "simple_network_example.csv",
        "name": "simple",
        "description": "Minimal: 11 edges with just 3 required columns"
    },
    {
        "file": "gene_regulatory_network_example.csv",
        "name": "genes",
        "description": "Gene regulation: 18 TF-gene interactions"
    }
]

for ex in examples:
    filepath = data_dir / ex["file"]
    if not filepath.exists():
        print(f"❌ File not found: {filepath}")
        continue
    
    df = pd.read_csv(filepath)
    project.add_datasource(ex["name"], dataframe=df, add_to_view=None)
    
    print(f"\n✅ Loaded: {ex['name']}")
    print(f"   File: {ex['file']}")
    print(f"   {ex['description']}")
    print(f"   Rows: {len(df)}")
    print(f"   Columns: {', '.join(df.columns)}")

# Create views with instructions for each dataset
project.set_view("1. Ligand-Receptor (Full)", {
    "initialCharts": {
        "ligand_receptor": [{
            "type": "text_box_chart",
            "param": [],
            "text": """
<h2>Ligand-Receptor Network (Full Example)</h2>

<p><b>This demonstrates ALL visual encodings!</b></p>

<h3>To visualize:</h3>
<ol>
  <li>Click <b>"Add Chart"</b> button</li>
  <li>Select datasource: <code>ligand_receptor</code></li>
  <li>Chart type: <b>"Network Graph (Flexible)"</b></li>
  <li>Select columns:
    <ul>
      <li><b>Source Node:</b> <code>source_cell</code></li>
      <li><b>Target Node:</b> <code>target_cell</code></li>
      <li><b>Edge Weight:</b> <code>interaction_score</code></li>
      <li><b>Category Filter:</b> <code>ligand_type</code></li>
      <li><b>Edge Length:</b> <code>spatial_distance</code></li>
      <li><b>Edge Color:</b> <code>pvalue</code></li>
      <li><b>Node Size:</b> <code>cell_count</code></li>
      <li><b>Node Type:</b> <code>cell_type</code></li>
    </ul>
  </li>
  <li>Click <b>"Add Chart"</b></li>
</ol>

<h3>What you'll see:</h3>
<ul>
  <li><b>Link thickness</b> = interaction strength (0.65-0.96)</li>
  <li><b>Link length</b> = spatial distance (8-31 µm)</li>
  <li><b>Link color</b> = p-value (red = significant)</li>
  <li><b>Node size</b> = cell count (9-35 interactions)</li>
  <li><b>Node color</b> = cell type (5 types)</li>
</ul>

<h3>The network shows:</h3>
<ul>
  <li>30 interactions between 20 cells</li>
  <li>7 ligand types: VEGFA, TGFB1, PDGFB, FGF2, TNF, IL6, IFNG</li>
  <li>5 cell types: T_cell, Cancer, Fibroblast, Endothelial, Macrophage</li>
  <li>Tumor microenvironment signaling pathways</li>
</ul>

<p><b>Try:</b> After creating, use filters to show only VEGFA or TGFB1 interactions!</p>
            """,
            "position": [5, 5],
            "size": [600, 650]
        }]
    }
})

project.set_view("2. Simple Network (Minimal)", {
    "initialCharts": {
        "simple": [{
            "type": "text_box_chart",
            "param": [],
            "text": """
<h2>Simple Network (Minimal Example)</h2>

<p><b>This shows the MINIMUM required columns!</b></p>

<h3>To visualize:</h3>
<ol>
  <li>Click <b>"Add Chart"</b> button</li>
  <li>Select datasource: <code>simple</code></li>
  <li>Chart type: <b>"Network Graph (Flexible)"</b></li>
  <li>Select columns:
    <ul>
      <li><b>Source Node:</b> <code>from_node</code></li>
      <li><b>Target Node:</b> <code>to_node</code></li>
      <li><b>Edge Weight:</b> <code>weight</code></li>
      <li><i>Leave all optional columns empty!</i></li>
    </ul>
  </li>
  <li>Click <b>"Add Chart"</b></li>
</ol>

<h3>What you'll see:</h3>
<ul>
  <li>11 edges connecting 7 nodes (A through G)</li>
  <li>Link thickness varies by weight (0.74-0.92)</li>
  <li>Clean, simple network layout</li>
  <li>No colors or sizes (uses defaults)</li>
</ul>

<p><b>This is the EASIEST way to visualize any network!</b></p>
<p>Just 3 columns: source, target, weight.</p>

<h3>Try it with your own data:</h3>
<p>Any CSV with 3 columns will work:
<ul>
  <li>Column 1: Source node (text)</li>
  <li>Column 2: Target node (text)</li>
  <li>Column 3: Edge weight (number)</li>
</ul>
</p>
            """,
            "position": [5, 5],
            "size": [600, 550]
        }]
    }
})

project.set_view("3. Gene Regulatory Network", {
    "initialCharts": {
        "genes": [{
            "type": "text_box_chart",
            "param": [],
            "text": """
<h2>Gene Regulatory Network</h2>

<p><b>Transcription factor → target gene interactions</b></p>

<h3>To visualize:</h3>
<ol>
  <li>Click <b>"Add Chart"</b> button</li>
  <li>Select datasource: <code>genes</code></li>
  <li>Chart type: <b>"Network Graph (Flexible)"</b></li>
  <li>Select columns:
    <ul>
      <li><b>Source Node:</b> <code>transcription_factor</code></li>
      <li><b>Target Node:</b> <code>target_gene</code></li>
      <li><b>Edge Weight:</b> <code>regulatory_score</code></li>
      <li><b>Category Filter:</b> <code>regulation_type</code> (optional)</li>
      <li><b>Edge Color:</b> <code>binding_affinity</code> (optional)</li>
    </ul>
  </li>
  <li>Click <b>"Add Chart"</b></li>
</ol>

<h3>What you'll see:</h3>
<ul>
  <li>18 regulatory interactions</li>
  <li>5 transcription factors: TP53, STAT3, NFKB, MYC, HIF1A</li>
  <li>Link thickness = regulatory strength (0.85-0.96)</li>
  <li>Link color = binding affinity (ChIP-seq scores)</li>
</ul>

<h3>After creating the chart:</h3>
<ol>
  <li>Click the ⚙️ <b>Settings</b> icon on the chart</li>
  <li>Enable <b>"Show Directionality"</b></li>
  <li>You'll see arrows: TF → gene</li>
</ol>

<p><b>Use case:</b> Visualize transcriptional regulatory networks, pathway analysis, gene circuits.</p>
            """,
            "position": [5, 5],
            "size": [600, 600]
        }]
    }
})

# Views save automatically when set_view is called

print("\n" + "="*70)
print("✅ PROJECT CREATED SUCCESSFULLY!")
print("="*70)
print(f"\nProject location: {project.dir}")
print("\nDatasets available:")
print("  1. ligand_receptor - 30 interactions, all visual features")
print("  2. simple - 11 edges, minimal example")
print("  3. genes - 18 TF-gene interactions")

print("\n" + "="*70)
print("NEXT STEPS - OPEN IN MDV GUI:")
print("="*70)
print("\n1. Launch MDV")
print("2. Load project: 'network_examples'")
print("3. You'll see 3 views - visit each one")
print("4. Follow the instructions in each view")
print("5. Click 'Add Chart' and select the columns as described")
print("6. Explore the interactive networks!")

print("\n" + "="*70)
print("QUICK TIPS:")
print("="*70)
print("\n📊 Start with view 2 (Simple Network)")
print("   - Easiest to understand")
print("   - Only 3 columns")
print("   - See the basics")

print("\n🧬 Then try view 3 (Gene Regulatory)")
print("   - Biological example")
print("   - Enable directionality arrows")
print("   - See TF → gene regulation")

print("\n🎯 Finally view 1 (Ligand-Receptor)")
print("   - Full-featured example")
print("   - All 8 columns used")
print("   - See every visual encoding")

print("\n💡 In each chart:")
print("   - Drag nodes to rearrange")
print("   - Click edges to highlight")
print("   - Use ⚙️ settings to adjust physics")
print("   - Filter data to explore subsets")

print("\n" + "="*70)
print("Ready to visualize! Open MDV now!")
print("="*70)

