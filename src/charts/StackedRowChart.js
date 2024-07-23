import {select} from "d3-selection";
import {easeLinear} from "d3-ease";
import SVGChart from "./SVGChart.js";
import BaseChart from "./BaseChart.js";
import { getHierarchicalNodes } from "../utilities/clustering.js";

class StackedRowChart extends SVGChart{
    constructor(dataStore,div,config){
        const axis_types={
            x:{type:"linear"},
            y:{type:"band"},
            ry:{}
        }
        if (!config.color_legend){
            config.color_legend={display:false};
        }
		super(dataStore,div,config,axis_types);
        this.config.type = "stacked_row_chart";
        this.dim = this.dataStore.getDimension("category_dimension");
        this.addToolTip();
        //redraw the chart
        this.onDataFiltered(null);
	}

    onDataFiltered(dim){
        //no need to change anything
        if (this.dim === dim || this.isPinned){
            return;
        }
        this.dim.getSankeyData(data=>{
            this.data = data;
            this.sortData();
            this.drawChart();            
        },this.config.param,{method:"stacked"})
    }

    setSize(x,y){
        super.setSize(x,y); 
        this.drawChart();
    }

    sortData(){
        this.sortedData= this.data.slice(0);
        const s = this.config.sort;
        if (s==="size" || s=== null){
            this.sortedData.sort((a,b)=>a.total-b.total)
        }
        else if (s==="name"){
            const v = this.dataStore.getColumnValues(this.config.param[0]);
            this.sortedData.sort((a,b)=>v[a.id].localeCompare(v[b.id]));
        }
        else if (s==="composition"){
            const cdata= this.data.map(x=>{
                const li = x.values.map(y=>y.per);
                li._id= x.id;
                return li;
            });
            const info = getHierarchicalNodes(cdata);
            this.sortedData= info.order.map(x=>this.data[x]);
        }   
    }

    getColorLegend(){
        return this.dataStore.getColorLegend(this.config.param[1]);
    }

    getSettings(){
        const settings=super.getSettings();
        const c = this.config;
        return settings.concat([
            {
                label:"Show Color Legend",
                type:"check",  
                current_value:c.color_legend?c.color_legend.display:true, 
                func:(x)=>{
                    c.color_legend.display=x;
                    this.setColorLegend();
                }
            },
            {
                type:"radiobuttons",
                label:"Sort Order",
                current_value:c.sort || "default",
                choices:[["Default","default"],["Size","size"],["Name","name"],["Composition","composition"]],             
                func:(v)=>{
                    c.sort=v;
                    this.sortData();
                    this.drawChart();      
                }
            }
        ]);
    }

    drawChart(tTime=500){
        const box = this._getContentDimensions();
        const c = this.config
        const vgap=4;
        const bnumber= this.sortedData.length;
        const bheight = (box.height-(vgap*(bnumber+1)))/bnumber;
        const colors= this.dataStore.getColumnColors(c.param[1]);
        const yvalues = this.dataStore.getColumnValues(c.param[0]);
        const xvalues = this.dataStore.getColumnValues(c.param[1]);
        
        const trans =  select(this.contentDiv).transition()
        .duration(tTime).ease(easeLinear);

        this.x_scale.domain([0,1]);
        this.y_scale.domain(this.sortedData.map(x=>yvalues[x.id]))
        this.updateAxis();

        this.graph_area.selectAll(".stacked-row")
        .data(this.sortedData, d => d.id)
        .join(
            enter => enter.append("g")
                .attr("transform",(d,i)=>`translate(0,${(vgap*(i+1))+(bheight*i)})`)
                .attr("class","stacked-row"),
            update=>update.transition(trans)
                .attr("transform",(d,i)=>`translate(0,${(vgap*(i+1))+(bheight*i)})`)
            
        )
        .selectAll("rect")
        .data(d=>{
            return d.values;
        },d=>d.id)
        .join(
            enter=>enter.append("rect")
                .attr("x",d=>d.perpos*box.width)
                .attr("y",0)
                .attr("width",d=>d.per*box.width)
                .attr("height",bheight)
                .attr("fill",d=>colors[d.id])
                .on("mouseover pointermove",(e,d)=>{           
                    this.showToolTip(e,`${xvalues[d.id]}<br>${d.count}`);
                }).
                on("mouseleave",()=>{
                    this.hideToolTip();
                }),
            update=>update.transition(trans)
                .attr("x",d=>d.perpos*box.width)
                .attr("width",d=>d.per*box.width)
                .attr("height",bheight)  
                .attr("fill",d=>colors[d.id])  //if color scheme changed           
        );
    }
}

BaseChart.types["stacked_row_chart"]={
    "class":StackedRowChart,
    name:"Stacked Row Chart",
    params:[
        {
            type:"text",
            name:"Category y axis"
        },
        {
            type:"text",
            name:"Category x axis"
        }
    ]
}

export default StackedRowChart;