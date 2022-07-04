import BaseChart from "./BaseChart.js";
import SVGChart from "./SVGChart.js";
import {select} from "d3-selection";
import { geoPath } from "d3-geo";
import {easeLinear} from "d3-ease";
import {scaleSequential,scaleLinear} from 'd3-scale';
import {interpolateYlGnBu} from 'd3-scale-chromatic';


class DensityPlot extends SVGChart{
    constructor(dataStore,div,config){
		const x_name= dataStore.getColumnName(config.param[0]);
		const y_name = dataStore.getColumnName(config.param[1]);
		if (!config.axis){

            config.axis={
                x:{size:30,label:x_name,textsize:13},
                y:{size:45,label:y_name,textsize:13}
            }
        }
		if (!config.title){
			config.title= x_name+" x "+y_name;	
		}

    	super(dataStore,div,config,{x:{},y:{}});
        const c = this.config;
        this.dim=this.dataStore.getDimension("density_dimension");
        const mm1= this.dataStore.getMinMaxForColumn(c.param[0]);
        const mm2= this.dataStore.getMinMaxForColumn(c.param[1]);
        this.contentDiv.addEventListener("wheel",e=>this.zoom(e));
        this.x_sc=1;
        this.y_sc=1;
        this._setRatio();
        this.x_scale.domain(this._getExtra(mm1));     
        this.y_scale.domain(this._getExtra(mm2,true));
        this.orig_y_scale=this.y_scale.domain();
        this.orig_x_scale= this.x_scale.domain();
        this.updateAxis();
        this.onDataFiltered();
    }

    _setRatio(){
        const p= this.config.param;
        const mm1= this.dataStore.getMinMaxForColumn(p[0]);
        const mm2= this.dataStore.getMinMaxForColumn(p[1]);
        const w =mm1[1]-mm1[0];
        const h =mm2[1]-mm2[0];
        this.whRatio=w/h;
    }

    _getExtra(mm,reverse=false){
        const margin = (mm[1]-mm[0])/10;
        if (reverse){
            return [mm[1]+margin,mm[0]-margin,];
        }
        return [mm[0]-margin,mm[1]+margin];
    }

    setSize(x,y){
        super.setSize(x.y);
        this.updateAxis();
        this.onDataFiltered()
    }

    zoom(e){
        var rect = this.contentDiv.getBoundingClientRect();
        const box =this._getContentDimensions()
    	const x= e.clientX-rect.left-box.left;
        const y = e.clientY-rect.top-box.top;
        const gx  = this.x_scale.invert(x);
        const gy = this.y_scale.invert(y);
        const factor = e.deltaY<0?0.9:1.1;
        const xd  = this.x_scale.domain();
        let new_len = (xd[1]-xd[0])*factor;
        let new_start = gx-((x/box.width)*new_len);
        this.x_scale.domain([new_start,new_start+new_len]);
        const yd= this.y_scale.domain();
        new_len = (yd[0]-yd[1])*factor;
        new_start = ((y/box.height)*new_len)+gy;
        this.y_scale.domain([new_start,new_start-new_len]);
        this._rescaleSVG();


       


        this.updateAxis();

       // this.onDataFiltered();
       // this.graph_area.attr("transform",`scale(${this.y_sc},${this.x_sc})`)


    }

    //rescales/offsets the g element conataining the chart based on the current axis
    _rescaleSVG(){
        const xr= this.x_scale.domain();
        const xlen = xr[1]-xr[0];
        const jj = this.x_scale.range()[1]/this.or_x_scale[1];
        const nxs  =(this.orig_x_scale[1]-this.orig_x_scale[0])/xlen*jj;
        const nsx = scaleLinear().domain(xr).range(this.or_x_scale);
        const nxp = nsx(xr[0])-nsx(this.orig_x_scale[0])*jj;
        const yr= this.y_scale.domain();
        const ylen = yr[0]-yr[1];
        const ll = this.y_scale.range()[1]/this.or_y_scale[1];
        const nys  =(this.orig_y_scale[0]-this.orig_y_scale[1])/ylen*ll;
        const nsy = scaleLinear().domain(yr).range(this.or_y_scale);
        const nyp = nsy(yr[0])-nsy(this.orig_y_scale[0])*ll;
        const x=this.margins.left;
        const y= this.margins.top;
        this.graph_area.attr("transform",`translate(${x-nxp},${y-nyp}) scale(${nxs},${nys})  `);
    }

 

    onDataFiltered(dim){
        if (this.dim === dim || this.isPinned){
            return;
        }
        if (dim === "all_removed"){
            this.filter=[];
            this.resetButton.style.display="none";
        }
   
     

        this.or_y_scale=[0,400/this.whRatio]; //this.y_scale.range();
        this.or_x_scale=[0,400]; //this.x_scale.range();



       
        const config={
            category: this.dataStore.getColumnValues(this.config.param[2]).indexOf(this.config.category),
            yscale:[this.orig_y_scale,this.or_y_scale],
            xscale:[this.orig_x_scale,this.or_x_scale]
            
        };

        this.margins.left;
        this.margins.right;
      
        this.dim.getDensityContours(data=>{
            this.data = data;
            this.drawChart();
            const x_sc= 1;//this.x_scale.range()[1]/x_range[1];
            const y_sc=1;//this.y_scale.range()[1]/y_range[1];
            const x=this.margins.left;
            const y= this.margins.top;
      

       
          //  this.graph_area.attr("transform",`translate(${x},${y}) scale(${1},${1})  `);
       
            
           // this.graph_area.attr("transform",`translate(${x},${y}) scale(${1},${1}) `);  
           this._rescaleSVG()

           
            
                   
         },this.config.param,config);
        
    }

    drawChart(tTime=400){

        var color = scaleSequential(interpolateYlGnBu)
    .domain([0, 0.5]); 
        const trans =  select(this.contentDiv).transition()
        .duration(tTime).ease(easeLinear);
     
        this.graph_area.selectAll(".dp-polygon")
        .data(this.data)
        .join("path")
        .attr("class","dp-polygon")
      .attr("d", geoPath())
      .attr("fill", function(d) { return color(d.value); })

    }


}

BaseChart.types["density_plot"]={
    name:"Density Plot",
    class:DensityPlot,
    params:[{
        type:"number",
        name:" x-axis"
    },
    {
        type:"number",
        name:"y axis"
    },
    {
        type:"text",
        name:"category"
    }
    ]
}