import SVGChart from "./SVGChart.js";
import BaseChart from "./BaseChart";
import {
    forceSimulation,
    forceLink,
    forceManyBody,
    forceCenter,
    drag,
    scaleSqrt,
    scaleLinear,
    schemeReds,
    zoom,
    zoomIdentity,
    pie,
    arc,
    select,
} from "d3";
import { getColorLegendCustom } from "../utilities/Color.js";
import { loadColumnData } from "@/datastore/decorateColumnMethod";

const color_schemes = {
    "blue yellow": [
        "#115f9a", "#1984c5", "#22a7f0", "#48b5c4",
        "#76c68f", "#a6d75b", "#c9e52f", "#d0ee11", "#d0f400",
    ],
    red: schemeReds[8].slice(1),
    "blue yellow red": ["blue", "yellow", "red"],
};

/**
 * FlexibleNetworkChart - Network visualization without metadata requirements
 * 
 * User selects all columns directly in the Add Chart dialog.
 * Works with ANY datasource that has interaction/network data.
 * 
 * param[0]: source_node (text) - Source node ID
 * param[1]: target_node (text) - Target node ID
 * param[2]: edge_weight (number) - Edge weight/score for thickness
 * param[3]: category_filter (text, optional) - Category to filter by
 * param[4]: edge_length (number, optional) - Value for link length
 * param[5]: edge_color (number, optional) - Value for link color
 * param[6]: node_size (number, optional) - Value for node size
 * param[7]: node_type (text, optional) - Node type for coloring
 */
class FlexibleNetworkChart extends SVGChart {
    constructor(dataStore, div, config) {
        config.link_strength = config.link_strength || 0.5;
        config.node_repulsion = config.node_repulsion || -400;
        config.title = config.title || "Network Graph";
        
        super(dataStore, div, config, {});
        
        const c = this.config;
        
        // Initialize scales
        this.linkThicknessScale = scaleLinear();
        this.nodeScale = scaleSqrt();
        
        // Set up link thickness (always from param[2])
        const weightQuantile = this.dataStore.getColumnQuantile(c.param[2], 0.01);
        this.linkThicknessScale.domain(weightQuantile).range([1, 10]);
        
        // Set up link length (param[4] if provided, else constant)
        if (c.param[4]) {
            this.linkLengthScale = scaleLinear();
            const lengthQuantile = this.dataStore.getColumnQuantile(c.param[4], 0.01);
            this.linkLengthScale.domain(lengthQuantile).range([30, 120]);
        }
        
        // Set up link color (param[5] if provided)
        if (c.param[5]) {
            this.linkColorScale = scaleLinear();
            const colorQuantile = this.dataStore.getColumnQuantile(c.param[5], 0.01);
            this.linkColorScale
                .domain([colorQuantile[0], (colorQuantile[0] + colorQuantile[1]) / 2, colorQuantile[1]])
                .range(["darkred", "orange", "lightgray"]);
        }
        
        // Set up node size (param[6] if provided, else constant)
        if (c.param[6]) {
            const sizeQuantile = this.dataStore.getMinMaxForColumn(c.param[6]);
            this.nodeScale.domain(sizeQuantile).range([5, 20]);
        } else {
            this.nodeScale.domain([1, 1]).range([8, 8]);
        }
        
        c.node_radius = c.node_radius || 8;
        c.show_labels = c.show_labels !== false;
        c.label_size = c.label_size || 10;
        c.show_directionality = c.show_directionality || false;
        c.link_opacity = c.link_opacity !== undefined ? c.link_opacity : 0.6;
        c.node_opacity = c.node_opacity !== undefined ? c.node_opacity : 1.0;
        c.use_pie_nodes = c.use_pie_nodes !== false; // Enable pie chart nodes by default
        
        // Initialize color legend
        if (!c.color_legend) {
            c.color_legend = { display: false };
        }
        
        // Set up node coloring
        this.nodeColorFunction = null;
        if (c.color_by) {
            const conf = {
                asArray: false,  // Return color strings for SVG, not RGB arrays
                overideValues: {
                    colorLogScale: c.log_color_scale,
                },
            };
            this._addTrimmedColor(c.color_by, conf);
            this.nodeColorFunction = this.dataStore.getColorFunction(c.color_by, conf);
            c.color_legend.display = true; // Show legend when color_by is set
        } else if (c.param[7]) {
            c.color_legend.display = true; // Show legend when node type is set
        }
        
        // Set up force simulation
        this.forceLink = forceLink().id((d) => d.id);
        this.forceManyBody = forceManyBody().distanceMax(200);
        this.simulation = forceSimulation()
            .force("link", this.forceLink.strength(c.link_strength))
            .force("charge", this.forceManyBody.strength(c.node_repulsion))
            .force("center", forceCenter(this.width / 2, this.height / 2));
        
        // Set up zoom behavior
        this.zoomBehavior = zoom()
            .scaleExtent([0.1, 10])  // Allow zoom from 10% to 1000%
            .on("zoom", (event) => {
                this.graph_area.attr("transform", event.transform);
            });
        
        this.svg.call(this.zoomBehavior);
        
        this.setColorLegend();
        this.reCalculate();
    }
    
    reCalculate() {
        const c = this.config;
        const index = this.dataStore.columnIndex;
        
        // Get required columns
        const sourceCol = index[c.param[0]];
        const targetCol = index[c.param[1]];
        const weightCol = index[c.param[2]];
        
        // Optional columns
        const categoryCol = c.param[3] ? index[c.param[3]] : null;
        const nodeSizeCol = c.param[6] ? index[c.param[6]] : null;
        const nodeTypeCol = c.param[7] ? index[c.param[7]] : null;
        
        // Build network
        const nodeSet = new Set();
        const nodeMetadata = {};
        const nodeComposition = {}; // Track type composition for each node (for pie charts)
        const nodeFilterStatus = {}; // Track if nodes have filtered connections
        this.linkData = [];
        const links = {}; // Track bidirectional links
        
        const f = this.dataStore.filterArray;
        
        // First pass: track all nodes and build composition data
        for (let i = 0; i < this.dataStore.size; i++) {
            const sourceId = sourceCol.values 
                ? sourceCol.values[sourceCol.data[i]]
                : sourceCol.data[i].toString();
            const targetId = targetCol.values
                ? targetCol.values[targetCol.data[i]]
                : targetCol.data[i].toString();
            
            if (sourceId === targetId) continue;
            
            // Track if this connection is filtered
            const isFiltered = f[i] > 0;
            
            // Mark nodes that have filtered connections
            if (isFiltered) {
                nodeFilterStatus[sourceId] = (nodeFilterStatus[sourceId] || 0) + 1;
                nodeFilterStatus[targetId] = (nodeFilterStatus[targetId] || 0) + 1;
            }
            
            // Build composition for pie charts (track all type occurrences)
            if (nodeTypeCol) {
                const nodeType = nodeTypeCol.values 
                    ? nodeTypeCol.values[nodeTypeCol.data[i]] 
                    : nodeTypeCol.data[i];
                
                if (!nodeComposition[sourceId]) {
                    nodeComposition[sourceId] = {};
                }
                nodeComposition[sourceId][nodeType] = (nodeComposition[sourceId][nodeType] || 0) + 1;
            }
        }
        
        // Second pass: build visible network
        for (let i = 0; i < this.dataStore.size; i++) {
            // Skip filtered rows
            if (f[i] > 0) continue;
            
            // Get source and target
            const sourceId = sourceCol.values 
                ? sourceCol.values[sourceCol.data[i]]
                : sourceCol.data[i].toString();
            const targetId = targetCol.values
                ? targetCol.values[targetCol.data[i]]
                : targetCol.data[i].toString();
            
            // Skip self-loops
            if (sourceId === targetId) continue;
            
            // Filter by category if specified
            if (c.category_filter && categoryCol) {
                const category = categoryCol.values 
                    ? categoryCol.values[categoryCol.data[i]]
                    : categoryCol.data[i].toString();
                if (category !== c.category_filter) continue;
            }
            
            // Handle bidirectional links - keep stronger direction
            const reverseKey = `${targetId}|${sourceId}`;
            if (links[reverseKey] !== undefined) {
                const reverseLink = this.linkData[links[reverseKey]];
                if (weightCol.data[i] > weightCol.data[reverseLink.d_index]) {
                    // Current is stronger - replace
                    reverseLink.source = sourceId;
                    reverseLink.target = targetId;
                    reverseLink.d_index = i;
                }
                continue;
            }
            
            nodeSet.add(sourceId);
            nodeSet.add(targetId);
            
            // Store node metadata
            if (!nodeMetadata[sourceId]) {
                nodeMetadata[sourceId] = {
                    size: nodeSizeCol ? nodeSizeCol.data[i] : 1,
                    type: nodeTypeCol ? (nodeTypeCol.values ? nodeTypeCol.values[nodeTypeCol.data[i]] : nodeTypeCol.data[i]) : null,
                    dataIndex: i  // Store data index for color_by
                };
            }
            if (!nodeMetadata[targetId]) {
                // For target nodes that don't appear as source, use the node ID as the type
                nodeMetadata[targetId] = {
                    size: 1,  // Default size for nodes that don't appear as source
                    type: nodeTypeCol ? targetId : null,  // Use node ID as its own type
                    dataIndex: undefined  // No direct data index for target-only nodes
                };
            }
            
            // Add link
            this.linkData.push({
                source: sourceId,
                target: targetId,
                d_index: i,
                weight: weightCol.data[i]
            });
            
            links[`${sourceId}|${targetId}`] = this.linkData.length - 1;
        }
        
        // Create node array with filter status and composition
        this.nodeData = Array.from(nodeSet).map(id => ({
            id,
            size: nodeMetadata[id]?.size || 1,
            type: nodeMetadata[id]?.type || null,
            dataIndex: nodeMetadata[id]?.dataIndex,
            hasFilteredConnections: (nodeFilterStatus[id] || 0) > 0,
            composition: nodeComposition[id] || null  // Composition data for pie charts
        }));
        
        this.drawChart();
    }
    
    drawChart() {
        // Clear only graph_area contents, not the entire SVG (preserves zoom)
        this.graph_area.selectAll("*").remove();
        const c = this.config;
        
        // Colors for nodes - prioritize color_by over node type
        let nodeColors = null;
        let useNumericColor = false;
        
        if (this.nodeColorFunction && c.color_by) {
            // Using numeric color_by column
            useNumericColor = true;
        } else if (c.param[7]) {
            // Using categorical node type - create value-to-color mapping
            const column = this.dataStore.columnIndex[c.param[7]];
            const colors = this.dataStore.getColumnColors(c.param[7]);
            nodeColors = {};
            column.values.forEach((value, index) => {
                nodeColors[value] = colors[index];
            });
        }
        
        // Create links
        const links = this.graph_area.append("g")
            .attr("class", "links")
            .selectAll("line")
            .data(this.linkData)
            .enter()
            .append("line")
            .attr("stroke-width", d => this.linkThicknessScale(d.weight))
            .style("stroke", d => {
                if (this.linkColorScale && c.param[5]) {
                    const val = this.dataStore.columnIndex[c.param[5]].data[d.d_index];
                    return this.linkColorScale(val);
                }
                return "currentcolor";
            })
            .style("opacity", c.link_opacity)
            .on("click", (e, d) => {
                this.dataStore.dataHighlighted([d.d_index], this);
            })
            .on("mouseover", (e, d) => {
                const weight = d.weight.toFixed(3);
                this.showToolTip(e, `${d.source.id} → ${d.target.id}<br>Weight: ${weight}`);
            })
            .on("mouseleave", () => {
                this.hideToolTip();
            });
        
        // Create arrows for directionality (if enabled)
        let arrows;
        if (c.show_directionality) {
            arrows = this.graph_area.append("g")
                .attr("class", "arrows")
                .selectAll("path")
                .data(this.linkData)
                .enter()
                .append("path")
                .attr("stroke-width", 2)
                .attr("fill", "none")
                .style("opacity", 1.0)
                .attr("stroke", d => {
                    if (this.linkColorScale && c.param[5]) {
                        const val = this.dataStore.columnIndex[c.param[5]].data[d.d_index];
                        return this.linkColorScale(val);
                    }
                    return "currentcolor";
                });
        }
        
        // Create nodes
        const nodes = this.graph_area.append("g")
            .attr("class", "nodes")
            .selectAll("g")
            .data(this.nodeData)
            .enter()
            .append("g");
        
        // Helper function to determine if a node should be a pie chart
        const shouldUsePie = (d) => {
            return c.use_pie_nodes && 
                   nodeColors && 
                   d.composition && 
                   Object.keys(d.composition).length > 1 &&
                   !useNumericColor; // Don't use pies for numeric color_by
        };
        
        // Create pie arc generator
        const pieGenerator = pie()
            .value(d => d.value)
            .sort(null);
        
        const arcGenerator = arc();
        
        // Add node visuals (either pie charts or simple circles)
        nodes.each((d, i, nodeElements) => {
            const nodeElement = nodeElements[i];
            const nodeRadius = c.param[6] ? this.nodeScale(d.size) : c.node_radius;
            
            if (shouldUsePie(d)) {
                // Draw pie chart
                arcGenerator
                    .innerRadius(0)
                    .outerRadius(nodeRadius);
                
                const pieData = Object.entries(d.composition).map(([type, count]) => ({
                    type: type,
                    value: count
                }));
                
                const arcs = pieGenerator(pieData);
                
                // Append pie slices
                const g = select(nodeElement);
                g.selectAll("path")
                    .data(arcs)
                    .enter()
                    .append("path")
                    .attr("d", arcGenerator)
                    .attr("fill", arc_d => nodeColors[arc_d.data.type] || "steelblue")
                    .attr("stroke", "#fff")
                    .attr("stroke-width", 1.5)
                    .style("fill-opacity", c.node_opacity);
            } else {
                // Draw simple circle
                select(nodeElement)
                    .append("circle")
                    .attr("r", nodeRadius)
                    .attr("fill", () => {
                        // Priority: color_by > node type > default
                        if (useNumericColor && d.dataIndex !== undefined) {
                            return this.nodeColorFunction(d.dataIndex);
                        }
                        if (nodeColors && d.type) {
                            return nodeColors[d.type] || "steelblue";
                        }
                        return "steelblue";
                    })
                    .attr("stroke", "#fff")
                    .attr("stroke-width", 1.5)
                    .style("fill-opacity", c.node_opacity);
            }
        });
        
        // Add click and drag handlers to all nodes
        nodes
            .on("click", (e, d) => {
                // Find all interactions for this node
                const indices = this.linkData
                    .filter(link => {
                        const sourceId = typeof link.source === 'object' ? link.source.id : link.source;
                        const targetId = typeof link.target === 'object' ? link.target.id : link.target;
                        return sourceId === d.id || targetId === d.id;
                    })
                    .map(link => link.d_index);
                this.dataStore.dataHighlighted(indices, this);
            })
            .call(drag()
                .on("start", (e, d) => this.dragstarted(e, d))
                .on("drag", (e, d) => this.dragged(e, d))
                .on("end", (e, d) => this.dragended(e, d)));
        
        // Add tooltips showing composition
        nodes.append("title")
            .text(d => {
                let tooltip = `${d.id}`;
                if (d.composition && Object.keys(d.composition).length > 0) {
                    tooltip += '\n\nComposition:';
                    const total = Object.values(d.composition).reduce((a, b) => a + b, 0);
                    Object.entries(d.composition)
                        .sort((a, b) => b[1] - a[1]) // Sort by count descending
                        .forEach(([type, count]) => {
                            const percent = ((count / total) * 100).toFixed(1);
                            tooltip += `\n${type}: ${count} (${percent}%)`;
                        });
                }
                if (d.hasFilteredConnections) {
                    tooltip += '\n\n(Has filtered connections)';
                }
                return tooltip;
            });
        
        // Add labels if enabled
        if (c.show_labels) {
            nodes.append("text")
                .attr("dx", d => {
                    // Position label outside node radius + small gap
                    const nodeRadius = c.param[6] ? this.nodeScale(d.size) : c.node_radius;
                    return nodeRadius + 4;
                })
                .attr("dy", 4)
                .style("font-size", `${c.label_size}px`)
                .style("fill", "currentcolor")
                .style("fill-opacity", c.node_opacity)
                .style("pointer-events", "none")
                .text(d => d.id);
        }
        
        // Set up simulation
        this.simulation.nodes(this.nodeData);
        
        if (this.linkLengthScale && c.param[4]) {
            const lengthCol = this.dataStore.columnIndex[c.param[4]];
            this.forceLink
                .links(this.linkData)
                .distance(d => this.linkLengthScale(lengthCol.data[d.d_index]));
        } else {
            this.forceLink.links(this.linkData).distance(50);
        }
        
        this.simulation.on("tick", () => {
            links
                .attr("x1", d => d.source.x)
                .attr("y1", d => d.source.y)
                .attr("x2", d => d.target.x)
                .attr("y2", d => d.target.y);
            
            if (arrows) {
                arrows.attr("d", (d) => {
                    let len = Math.sqrt(
                        (d.source.y - d.target.y) ** 2 + (d.target.x - d.source.x) ** 2,
                    );
                    len = len === 0 ? 0.1 : len;
                    const f = (len / 2 - 5) / len;
                    const x2 = (d.target.x - d.source.x) / 2 + d.source.x;
                    const y2 = (d.target.y - d.source.y) / 2 + d.source.y;
                    const xmido = (d.target.x - d.source.x) * f + d.source.x;
                    const ymido = (d.target.y - d.source.y) * f + d.source.y;
                    const xo = (5 / len) * (d.source.y - d.target.y);
                    const x1 = xmido + xo;
                    const x3 = xmido - xo;
                    const yo = (5 / len) * (d.target.x - d.source.x);
                    const y1 = ymido + yo;
                    const y3 = ymido - yo;
                    return `M ${x1} ${y1} ${x2} ${y2} ${x3} ${y3}`;
                });
            }
            
            nodes.attr("transform", d => `translate(${d.x},${d.y})`);
        });
        
        this.simulation.alpha(1).restart();
    }
    
    dragstarted(event, d) {
        if (!event.active) this.simulation.alphaTarget(0.3).restart();
        d.fx = d.x;
        d.fy = d.y;
    }
    
    dragged(event, d) {
        d.fx = event.x;
        d.fy = event.y;
    }
    
    dragended(event, d) {
        if (!event.active) this.simulation.alphaTarget(0);
        d.fx = null;
        d.fy = null;
    }
    
    @loadColumnData
    setParams(params) {
        this.config.param = params;
        const c = this.config;
        
        // Rebuild all scales with new columns
        // Set up link thickness (always from param[2])
        const weightQuantile = this.dataStore.getColumnQuantile(c.param[2], 0.01);
        this.linkThicknessScale.domain(weightQuantile).range([1, 10]);
        
        // Set up link length (param[4] if provided, else constant)
        if (c.param[4]) {
            if (!this.linkLengthScale) {
                this.linkLengthScale = scaleLinear();
            }
            const lengthQuantile = this.dataStore.getColumnQuantile(c.param[4], 0.01);
            this.linkLengthScale.domain(lengthQuantile).range([30, 120]);
        } else {
            this.linkLengthScale = null;
        }
        
        // Set up link color (param[5] if provided)
        if (c.param[5]) {
            if (!this.linkColorScale) {
                this.linkColorScale = scaleLinear();
            }
            const colorQuantile = this.dataStore.getColumnQuantile(c.param[5], 0.01);
            this.linkColorScale
                .domain([colorQuantile[0], (colorQuantile[0] + colorQuantile[1]) / 2, colorQuantile[1]])
                .range(["darkred", "orange", "lightgray"]);
        } else {
            this.linkColorScale = null;
        }
        
        // Set up node size (param[6] if provided, else constant)
        if (c.param[6]) {
            const sizeQuantile = this.dataStore.getMinMaxForColumn(c.param[6]);
            this.nodeScale.domain(sizeQuantile).range([5, 20]);
        } else {
            this.nodeScale.domain([1, 1]).range([8, 8]);
        }
        
        // Update color legend if node type changed
        if (c.param[7]) {
            c.color_legend.display = true;
        }
        
        this.setColorLegend();
        this.reCalculate();
    }
    
    onDataFiltered(dim) {
        this.reCalculate();
    }
    
    getColorLegend() {
        const c = this.config;
        // Priority: color_by > node type
        if (c.color_by) {
            return this.dataStore.getColorLegend(c.color_by, {
                name: this.dataStore.getColumnName(c.color_by)
            });
        }
        // Fall back to node type
        if (c.param[7]) {
            return this.dataStore.getColorLegend(c.param[7], {
                name: this.dataStore.getColumnName(c.param[7])
            });
        }
        return null;
    }
    
    colorByColumn(column) {
        this.config.color_by = column;
        const conf = {
            asArray: false,  // Return color strings for SVG, not RGB arrays
            overideValues: {
                colorLogScale: this.config.log_color_scale,
            },
        };
        this._addTrimmedColor(column, conf);
        this.nodeColorFunction = this.dataStore.getColorFunction(column, conf);
        
        if (!this.config.color_legend) {
            this.config.color_legend = { display: true };
        }
        this.config.color_legend.display = true;
        
        this.setColorLegend();
        this.drawChart();
    }
    
    colorByDefault() {
        this.config.color_by = undefined;
        this.nodeColorFunction = null;
        
        // Keep legend if node type is selected
        if (!this.config.param[7]) {
            if (this.legend) {
                this.legend.remove();
                this.legend = undefined;
            }
        } else {
            this.setColorLegend();
        }
        
        this.drawChart();
    }
    
    getColorOptions() {
        return {
            colorby: ["integer", "double", "int32"],  // Support numeric columns for coloring
            // Note: color_overlay removed - not needed for network charts and causes unnecessary redraws
        };
    }
    
    remove() {
        this.simulation.on("tick", null);
        super.remove();
    }
    
    setSize(x, y) {
        super.setSize(x, y);
        this.simulation.force(
            "center",
            forceCenter(this.width / 2, this.height / 2),
        );
        this.simulation.alpha(0.3).restart();
    }
    
    resetZoom() {
        this.svg.transition()
            .duration(750)
            .call(this.zoomBehavior.transform, zoomIdentity);
    }
    
    getSettings() {
        const settings = super.getSettings();
        const c = this.config;
        
        settings.push({
            type: "slider",
            max: 2,
            min: 0,
            current_value: c.link_strength,
            label: "Link Strength",
            func: (x) => {
                c.link_strength = x;
                this.forceLink.strength(x);
                this.simulation.alpha(0.3).restart();
            },
        });
        
        settings.push({
            type: "slider",
            max: 0,
            min: -1000,
            current_value: c.node_repulsion,
            label: "Node Repulsion",
            func: (x) => {
                c.node_repulsion = x;
                this.forceManyBody.strength(x);
                this.simulation.alpha(0.3).restart();
            },
        });
        
        settings.push({
            type: "slider",
            max: 30,
            min: 3,
            current_value: c.node_radius,
            label: "Node Radius",
            func: (x) => {
                c.node_radius = x;
                if (!c.param[6]) {
                    // No node size column - update fixed radius scale
                    this.nodeScale.domain([1, 1]).range([x, x]);
                } else {
                    // Has node size column - adjust the scale range
                    const domain = this.nodeScale.domain();
                    this.nodeScale.domain(domain).range([x * 0.5, x * 2]);
                }
                // Redraw chart to update both circles and pie chart nodes
                this.drawChart();
            },
        });
        
        settings.push({
            type: "check",
            current_value: c.show_labels,
            label: "Show Node Labels",
            func: (x) => {
                c.show_labels = x;
                this.drawChart();
            },
        });
        
        settings.push({
            type: "slider",
            max: 20,
            min: 5,
            current_value: c.label_size,
            label: "Label Size",
            func: (x) => {
                c.label_size = x;
                this.svg.selectAll("text")
                    .style("font-size", `${x}px`);
            },
        });
        
        // Only show pie chart option if node type column is selected
        if (c.param[7] && !c.color_by) {
            settings.push({
                type: "check",
                current_value: c.use_pie_nodes,
                label: "Show Composition as Pie Charts",
                func: (x) => {
                    c.use_pie_nodes = x;
                    this.drawChart();
                },
            });
        }
        
        settings.push({
            type: "slider",
            max: 1,
            min: 0.1,
            current_value: c.link_opacity,
            label: "Link Opacity",
            func: (x) => {
                c.link_opacity = x;
                this.svg.selectAll(".links line")
                    .style("opacity", x);
            },
        });
        
        settings.push({
            type: "slider",
            max: 1,
            min: 0.1,
            current_value: c.node_opacity,
            label: "Node Opacity",
            func: (x) => {
                c.node_opacity = x;
                // Update both circles and pie chart paths
                this.svg.selectAll(".nodes circle").style("fill-opacity", x);
                this.svg.selectAll(".nodes path").style("fill-opacity", x);
            },
        });
        
        settings.push({
            type: "check",
            current_value: c.show_directionality,
            label: "Show Directionality",
            func: (x) => {
                c.show_directionality = x;
                this.drawChart();
            },
        });
        
        // Category filter dropdown (only if category column is selected)
        if (c.param[3]) {
            const categoryValues = this.dataStore.getColumnValues(c.param[3]);
            const categories = categoryValues.map((val) => {
                return { name: val, value: val };
            });
            // Add "All" option at the beginning
            categories.unshift({ name: "All", value: null });
            
            settings.push({
                type: "dropdown",
                label: `Filter by ${this.dataStore.getColumnName(c.param[3])}`,
                values: [categories, "name", "value"],
                current_value: c.category_filter || null,
                func: (x) => {
                    c.category_filter = x;
                    this.reCalculate();
                },
            });
        }
        
        // Show legend toggle only if node type column is selected AND not using color_by
        if (c.param[7] && !c.color_by) {
            settings.push({
                type: "check",
                current_value: c.color_legend.display,
                label: "Show Node Type Legend",
                func: (x) => {
                    c.color_legend.display = x;
                    this.setColorLegend();
                },
            });
        }
        
        // Reset zoom button
        settings.push({
            type: "button",
            label: "Reset Zoom",
            func: () => {
                this.resetZoom();
            },
        });
        
        return settings;
    }
}

// Register as a new chart type - works with ANY datasource!
BaseChart.types["flexible_network_chart"] = {
    name: "Network Graph",
    class: FlexibleNetworkChart,
    
    // NO required metadata! Works with any table
    
    // Define parameters for GUI column selection
    params: [
        { type: "text", name: "Source Node" },
        { type: "text", name: "Target Node" },
        { type: "number", name: "Edge Weight" },
        { type: "text", name: "Category Filter (optional)" },
        { type: "number", name: "Edge Length (optional)" },
        { type: "number", name: "Edge Color (optional)" },
        { type: "number", name: "Node Size (optional)" },
        { type: "text", name: "Node Type (optional)" },
    ],
};

export default FlexibleNetworkChart;

