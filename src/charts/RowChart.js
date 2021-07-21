import {select} from "d3-selection";
import {easeLinear} from "d3-ease";
import {CategoryChart} from "./CategoryChart.js";
import { BaseChart } from "./BaseChart.js";

class RowChart extends CategoryChart{
    constructor(dataStore,div,config){
		super(dataStore,div,config,{x:{}});
        this.config.type = "row_chart";
        //redraw the chart
        this.onDataFiltered(null);
	}

    drawChart(tTime=400){
        const trans =  select(this.contentDiv).transition()
        .duration(tTime).ease(easeLinear);
        const chartWidth= this.width-this.margins.left-this.margins.right;
        const chartHeight = this.height-this.margins.bottom- this.margins.top;
        const nBars = this.rowData.length;
        const vals = this.dataStore.getColumnValues(this.config.param);
        const barHeight= (chartHeight-(nBars+1)*3)/nBars;
        
  
        this.x_scale.domain([0,this.maxCount]);
        this.updateAxis();
        const units = chartWidth/this.maxCount;
        const self = this;
        this.graph_area.selectAll(".row-bar")
        .data(this.rowData,(d,i)=>i)
        .join("rect")
        .on("click",(e,d)=>{
            self.filterCategories(vals[d[1]],e.shiftKey);
        })
        .transition(trans)
        .style("fill",(d,i)=>{
            if (self.filter.length>0){
                if (self.filter.indexOf(vals[i])===-1){
                    return "lightgray"
                }
            }
            return self.colors[i];
          
        })
        .attr("class","row-bar")
        .attr("x",0)
        .attr("width",d=>d[0]*units)
        .attr("y",(d,i)=>((i+1)*3)+(i*barHeight))
        .attr("height",barHeight)
    

        this.graph_area.selectAll(".row-text")
        
        .data(this.rowData,(d,i)=>i)
        .join("text")
        .on("click",(e,d)=>{
            self.filterCategories(vals[d[1]],e.shiftKey);
        })
        .attr("class","row-text")
        .transition(trans)
        .text((d,i)=>vals[i])
        .attr("x",5)
        .attr("y",(d,i)=>((i+1)*3)+(i*barHeight)+(barHeight/2))
        //.attr("text-anchor", "middle")
        .attr("dominant-baseline", "central") 

    }
}

BaseChart.types["row_chart"]={
    "class":RowChart,
    name:"Row Chart",
    params:[{
        type:"text",
        name:"Category"
    }]

}

export {RowChart};
