
import {select} from "d3-selection";
import {easeLinear} from "d3-ease";
import SVGChart from "./SVGChart.js";
import BaseChart from "./BaseChart.js";
import {brushX} from "d3-brush";

class HistogramChart extends SVGChart{
    constructor(dataStore,div,config){
        if (Array.isArray(config.param)){
            config.param=config.param[0];
        }
        const pname= dataStore.getColumnName(config.param);
		if (!config.axis){
            config.axis={
                x:{size:30,label:pname,textsize:13},
                y:{size:45,label:"frequency",textsize:13}
            }
        }
		if (!config.title){
			config.title=pname;	
		}
       
		super(dataStore,div,config,{x:{},y:{}});
        const c  =this.config
        c.type="bar_chart";
        
        
        		
        const mm = this.dataStore.getMinMaxForColumn(config.param);
        this._setConfigValue(c,"display_max",mm[1]);
        this._setConfigValue(c,"display_min",mm[0]);
        this._setConfigValue(c,"bin_number",10)
 
        this.brush = brushX()
            .on("end",(e)=>{
                const s = e.selection;
                if (s==null){
                    return;
                }
                const offset = this.margins.left;
                this.setRange(this.x_scale.invert(s[0]-offset),this.x_scale.invert(s[1]-offset));            
            })


        this._setBrushExtent();
        this.dim = this.dataStore.getDimension("range_dimension");
        this.filter=null;
        this.binFiltered=new Array(this.config.bin_number).fill(false);
        this.onDataFiltered(null,true);
  
     
	}


    _setBrushExtent(){
        this.brush.extent([[this.margins.left,0],[this.width-this.margins.right,this.height-this.margins.bottom]]);
        this.svg.call(this.brush);
    }

    setSize(x,y){
        super.setSize(x,y);   
        this._setBrushExtent();
        this.drawChart();
        this.svg.call(this.brush.move,null);
       
             
    }

    remove(notify=true){
        this.dim.destroy(notify);
        super.remove();
    }

    


    removeFilter(){
        this.dim.removeFilter();
        this.filter=null;
        this.binFiltered=new Array(this.config.bin_number).fill(false);
        this.drawChart();
        this.svg.call(this.brush.move,null);

    }

    _calculateFilterBins(){
        const c = this.config;
        const nBins = c.bin_number;
        this.binFiltered=new Array(nBins).fill(false);
        if (!this.filter){
            return;
        }
        const iSize = (c.display_max-c.display_min)/(c.bin_number);
        let st  = c.display_min;
        for (let i=0;i<c.bin_number;i++){
            const pos=st+(i*iSize)
            if (pos<= this.filter[0] || pos>this.filter[1]){
                this.binFiltered[i] =true;
            }
        }    
    }

    setRange(start,end){
        this.filter=[start,end];      
        this._calculateFilterBins();
        this.dim.filter("filterRange",[this.config.param],{min:start,max:end});
        this.drawChart();
        this.resetButton.style.display = "inline";
    }

    getFilter(){
        if (!this.filter){
            return null;
        }
        const f={};
        f[this.config.param]=this.filter.slice(0);
        return f;
    }

   
    onDataFiltered(dim){
        //no need to change anything
        
        if (this.dim === dim || this.isPinned){
            return;
        }
        const bn = this.config.bin_number;
        const config=
        {             
            bins:this.config.bin_number,
            max:this.config.display_max,
            min:this.config.display_min
        }
        if (dim==="all_removed"){
            this.filter=null;
            this.binFiltered.fill(false);
            this.svg.call(this.brush.move,null);
            this.resetButton.style.display="none";
        }
       
    
        
        this.dim.getBins(data=>{
            this.originalBinData=data;
            this.binData=data.slice(0,bn);
            this.binData[bn-1]+=data[bn];
            this.drawChart();
        },this.config.param,config
       )
       
    }

    drawChart(){
        const trans = select(this.contentDiv).transition()
        .duration(400).ease(easeLinear);
        const cdim = this._getContentDimensions();
        const nBins = this.config.bin_number;
        const barWidth= (cdim.width-(nBins+1)*3)/nBins;
        let maxBar= Math.max(...this.binData);
        maxBar=maxBar?maxBar:1;

        this.x_scale.domain([this.config.display_min,this.config.display_max]);
		this.y_scale.domain([maxBar,0]);  

        this.updateAxis();
        const units = cdim.height/maxBar;
        const self = this;
        let join =this.graph_area.selectAll(".bar")
        
        .data(this.binData)
        .join("rect")       
        .transition(trans)
        .style("fill",(d,i)=>{
           return self.binFiltered[i]?"gray":"blue"
        })
        .attr("class","bar")
        .attr("x",(d,i)=>((i+1)*3)+(i*barWidth))
        .attr("width",barWidth)
        .attr("y",d=>cdim.height-(d*units))
        .attr("height",d=>d*units)

    }

    pinChart(){
        this.isPinned=true;
    }


    unpinChart(){
        this.isPinned = false;
        this.onDataFiltered();
    }

    getSettings(){
        const settings = super.getSettings();
        const c = this.config;
        const maxMin = this.dataStore.getMinMaxForColumn(c.param)
        return settings.concat([
          
            {
                type:"doubleslider",
                min:maxMin[0],
                max:maxMin[1],
                doc:this.__doc__ || document,
                current_value:[c.display_min,c.display_max],
                label:"max min values",
                info:{text:`Sometimes outliers can skew the histogram to one side.
                     Setting a min/max value will add all counts less/greater than the 
                     value specified to the first/last bin.`,position:"left"},
                func:(x,y)=>{
                    c.display_min=x;
                    c.display_max=y;
                    this._calculateFilterBins();
                    this.onDataFiltered();
                }
            },
            {
                type:"spinner",
                max:100,
                min:2,
                current_value:c.bin_number,
                label:"Number of bins",
                func:(x)=>{
                    c.bin_number=x;
                    this._calculateFilterBins();
                    this.onDataFiltered();
                }
            }
        ]);      
    }
}

BaseChart.types["bar_chart"]={
    name:"Histogram",
    class:HistogramChart,
    params:[{
        type:"number",
        name:"Frequency Data"
    }]
}

export default HistogramChart;