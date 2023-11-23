import SVGChart from "./SVGChart.js";
import {line,area} from "d3-shape";
import {easeLinear} from "d3-ease";
import {select} from "d3-selection";
import {scaleBand,scaleLinear} from "d3-scale";
import {axisLeft} from "d3-axis";
import BaseChart from "./BaseChart";
import { ClipSpace } from "@luma.gl/core";


class MultiLineChart extends SVGChart{
    constructor(dataStore,div,config){
		const x_name= dataStore.getColumnName(config.param[0]);
		if (!config.axis){

            config.axis={
                x:{size:30,label:x_name,textsize:13},
                y:{size:45,label:"density",textsize:13}
            }
        }
		if (!config.title){
			config.title= x_name;	
		}
        if (config.scaletrim && !(config.band_width)){
            config.band_width=0.1;
        }
        
    	super(dataStore,div,config,{x:{},y:{type:config.stacked?"band":"linear"}});
        this.config.type="multi_line_chart";

        const c = this.config;
        if (!c.color_legend ){ 
            c.color_legend={display:c.stacked?false:true}
        }
        //set default band width
        let mm = this.dataStore.getMinMaxForColumn(c.param[0]);
        let xx= mm.slice(0);
        let yy= [1,0];
        if (config.scaletrim){
            xx = [0,1];
        }
        this.x_scale.domain(xx)
        //this.y_scale.domain(yy);
        this.defaultBandWidth= (mm[1]-mm[0])/100;
        c.band_width =c.band_width || this.defaultBandWidth;
        c.intervals = c.intervals || 40;

        this.dim= this.dataStore.getDimension("catrange_dimension");
        this.data=[];
        this.addToolTip();
      
        this.setColorLegend();
        this.stackLines(c.stacked);
        this.scaleTrim(c.scaletrim);
        this.onDataFiltered(); 
    }

    stackLines(st){
        
        this.config.stacked=st;
        if (st){
            this.y_scale=scaleBand();
            this.y_axis_call= axisLeft(this.y_scale);

            this.config.color_legend.display=false;
            this.setColorLegend();
            this.y_scale.domain(this.dataStore.getColumnValues(this.config.param[1]));

        }
        else{
            this.y_scale=scaleLinear();
            this.y_axis_call= axisLeft(this.y_scale);
            this.y_scale.domain([1,0]);
        }
    }

    changeCategories(column,val1,val2){
        this.config.param[1]=column;
        this.onDataFiltered();
    }

    scaleTrim(val){
        const c = this.config;
        if (!(val) || val==="none"){
            delete c.scaletrim;
            const mm = this.dataStore.getMinMaxForColumn(c.param[0]);
            this.x_scale.domain(mm)
        }
        else{
            c.scaletrim=val+"";
            this.x_scale.domain([0,1])
            c.band_width=0.1;
        }

    }


    getColorLegend(){
        return this.dataStore.getColorLegend(this.config.param[1]);
    }

    onDataFiltered(dim){
        if (this.isPinned){
            return;
        }     
        else if (dim !== this.dim){
            if (dim === "all_removed"){
          
                this.resetButton.style.display="none";  
            }
            const c = this.config;
            this.ticks=this.x_scale.ticks(c.intervals);
            this.dim.getKernalDensity(data=>{
                this.data=data;
                this.drawChart();
            },[c.param[1],c.param[0]],
            {
                ticks:this.ticks,
                bandwidth:c.band_width,
                intervals:c.intervals,
                scaletrim:c.scaletrim
            })
        }       	
	}

    remove(notify=true){
        this.dim.destroy(notify);
        super.remove();
    }

   
    drawChart(){
        
        const c = this.config;
        const vals = this.dataStore.getColumnValues(c.param[1])
        const trans = select(this.contentDiv).transition()
        .duration(400).ease(easeLinear);
        const cdim = this._getContentDimensions();
        const y_scales=[];
        const display_data=[];
        const rowHeight = cdim.height/this.data.length;
        const axisd=[];
        for (let i=0;i<this.data.length;i++){
            const scale= c.stacked?rowHeight/this.data[i].max:(cdim.height/this.data[i].max);
            const arr= c.stacked?this.data[i].map(x=>((rowHeight*i)+rowHeight)-(x*scale)) :this.data[i].map(x=>cdim.height-(x*scale));
            arr.id= this.data[i].id;
            axisd.push(vals[arr.id])
            display_data.push(arr);
            
        }
        this.y_scale.domain(axisd);
        this.updateAxis();
        const xpos= this.ticks.map(x=>this.x_scale(x));
        const self = this;
        const colors = this.dataStore.getColumnColors(this.config.param[1]);
        const ga =this.graph_area.selectAll(".chart-line")
        .data(display_data);
        ga.join("path")     
        .attr("class","chart-line")
        .on("mouseover pointermove",(e,d)=>{          
            const label = vals[d.id];
            self.showToolTip(e,label);
        }).
        on("mouseout",()=>{
            self.hideToolTip();
        })
        .transition(trans)
        .attr("fill", "none")
        //.attr("opacity",0.7)

        .attr("stroke", (d,i)=>{
            
            return colors[d.id]
        })
        .attr("stroke-width", 2)
        .attr("stroke-linejoin", "round")
       
       
        .attr("d",line()
          .x((d,i)=>{
              return xpos[i]
          })
          .y((d,i)=>{
              return d;
            })
           
        )
     

      if (c.fill){
            const a= area()
            .x(function(d,i) { return xpos[i]; })
            .y0(function(d,i) { return d; })
            .y1((d,i)=>{
                c.stacked?(rowHeight*d.id)+rowHeight:cdim.height
            });
        
        this.graph_area.selectAll(".chart-area")
        .data(display_data)
        .join("path")     
        .attr("class","chart-area")
        .transition(trans)
        .attr("fill", d=>colors[d.id])
        .style("opacity", 0.5)
        .attr("d", function(d,i){
            return area()
            .x(function(dd,ii) { return xpos[ii]; })
            .y0(function(dd) { return dd; })
            .y1(d=>{
                return c.stacked?(rowHeight*i)+rowHeight:cdim.height;
            })(d);
        });
    
    }else{
        this.graph_area.selectAll(".chart-area").remove();
    }
  
    }

    unpinChart(){
        super.unpinChart();
        this.onDataFiltered();
    }



    setSize(x,y){
        super.setSize(x,y);
        this.drawChart();
    }

    
    getSettings(){
        const s = super.getSettings();
        const c = this.config;
        const mm = this.dataStore.getMinMaxForColumn(c.param[0]);
        const cols = this.dataStore.getColumnList("text");
       
        s.splice(1,0,
            {
                type:"dropdown",
                current_value:c.param[1],
                label:"Change Cetegories (y axis)",
                values:[cols,"name","field"],
                func:x=>{
                    this.changeCategories(x,"val1","val2");
                }


            },
            {          
                type:"check",
                current_value:c.stacked,
                label:"Stack Lines",
                func:(x)=>{             
                    this.stackLines(x);
                    this.drawChart();
                    
                }
            },
            {

            type:"check",
            current_value:c.fill,
            label:"Fill Lines",
            func:(x)=>{             
                c.fill=x;
                this.drawChart();
                
            }
        },
            {
                type:"radiobuttons",
                label:"Scale/Trim to Percentile",
                current_value:c.scaletrim || "none",
                choices:[["No Trim","none"],["0.001","0.001"],["0.01","0.01"],["0.05","0.05"]],             
                func:(v)=>{
                 
                    this.scaleTrim(v)
                    this.onDataFiltered();
                           
                }
            },
            {
            
                type:"slider",
                max:mm[1]/10,
                min:mm[0]<0?0:mm[0],
                doc:this.__doc__,
                current_value:c.band_width,
                label:"Band Width",
                func:(x)=>{
                    c.band_width=x;
                    this.onDataFiltered();
                }
            },
            {
            
                type:"slider",
                max:100,
                min:10,
                step:1,
                doc:this.__doc__,
                current_value:c.intervals,
                label:"Intervals",
                func:(x)=>{
                    c.intervals=x;
                    this.onDataFiltered();
                }
            },
            {
                label:"Show Color Legend",
                type:"check",
                current_value:c.color_legend?c.color_legend.display:true, 
                func:(x)=>{
                    c.color_legend.display=x;     
                    this.setColorLegend();
                }

            }
        );
        return s;
    }

}

BaseChart.types["multi_line_chart"]={
    "class":MultiLineChart,
    name:"Multi Line Chart",
    methodsUsingColumns:["changeCategories"],
    params:[
     
        {
            type:"number",
            name:"Value (X axis)"
        },
        {
            type:"text",
            name:"Categories to show"
        },
    ]

}

export default MultiLineChart;