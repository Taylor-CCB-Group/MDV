
import BaseChart from "./BaseChart.js"
import SVGChart from "./SVGChart.js";
import { parseNewick } from "../utilities/clustering.js";
import {easeLinear} from "d3-ease";
import {cluster} from "d3-hierarchy";
import {select} from "d3-selection";
import {max} from "d3"


class TreeDiagram extends SVGChart{
    constructor(dataStore,div,config){
		super(dataStore,div,config,{});
        const hld= this.dataStore.getHighlightedData();
        if (hld){
            this.onDataHighlighted({data:this.dataStore.rowData.get(hld[0])});
           
        }    
	}

    setBranchLength(d, y0, k) {
        d.y = (y0 += d.data.length) * k;
        if (d.children) d.children.forEach(d => this.setBranchLength(d, y0, k));
      }

   
    drawChart(){
        if (!this.treeData){
            return
        }
        const trans =  select(this.contentDiv).transition()
        .duration(500).ease(easeLinear);
        const b= this._getContentDimensions();
        const treeMap = cluster().size([b.height,b.width-50]).separation((a,b)=>1);
        
        const nodes= treeMap(this.treeData);
        function maxLength(d) {
            return d.data.length + (d.children ? max(d.children, maxLength) : 0);
        }
    
        this.setBranchLength(this.treeData, this.treeData.data.length = 0,(b.width-50)/maxLength(this.treeData));
        
        const desc=  nodes.descendants();
        const nml = max(desc,d=>d.data.length);
        this.graph_area.selectAll(".tree-link")
            .data(desc.slice(1))
            .join("path")
            .attr("class", "tree-link")
            .style("stroke","currentColor")
            .style("fill","none")
            .style("stroke-width","1px")
            .transition(trans)
            .attr("d", d=>`M ${d.y} ${d.x} L${d.parent.y} ${d.x} L${d.parent.y} ${d.parent.x}`);
        this.graph_area
            .selectAll("text")
            .data(this.treeData.leaves())
            .join("text")
            .text(d => d.data.name.replace(/_/g, " "))
            .attr("font-size","9px")
            .attr("alignment-baseline","middle")
            .style("stroke","currentColor")
            .transition(trans)
            .attr("x", d=>d.y+2)
            .attr("y", d=>d.x);         
    }

   
    setSize(x,y){
        super.setSize(x,y);
        this.drawChart();
    }

    onDataHighlighted(data){
       const nwk = data.data?.nwk_data;
       if (nwk){
            this.treeData= parseNewick(nwk);
            this.drawChart();
       }
    }

  
    getSettings(){
        return super.getSettings();     
    }
    
}

BaseChart.types["row_summary_box"]={
    "class":TreeDiagram,
    name:"Tree Diagram",
    required:["tree_diagram"],
    params:[]
}

export default TreeDiagram;