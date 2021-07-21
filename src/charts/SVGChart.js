import {BaseChart} from "./BaseChart.js";
import {select} from "d3-selection";
import {scaleLinear} from "d3-scale";
import "d3-transition";
import {axisLeft,axisBottom} from "d3-axis";


class SVGChart extends BaseChart{
    constructor(dataStore,div,config,axisTypes){
		super(dataStore,div,config);
        this.svg = select(this.contentDiv)
        .append("svg")
        .attr("width", this.width)
        .attr("height",this.height)
     
            
        this.config.axis =this.config.axis || {};
        const cax = this.config.axis;
       
        for (let at of [["x",25],["y",40],["ry",25],["tx",25]]){
            if (axisTypes[at[0]]){
                cax[at[0]] = cax[at[0]] || {textsize:13,label:"",size:at[1]}
            }

        }

        this._setMargins();
          

        this.graph_area=this.svg
            .append("g")
            .attr("transform", `translate(${this.margins.left},${this.margins.top})`);
        

      // add the x Axis
        if (cax.x){
            this.x_axis_svg=this.graph_area.append("g");
            this.x_scale = scaleLinear();
            this.x_axis_call= axisBottom(this.x_scale);
            this.x_axis_label=this.svg.append("text")
                .style("text-anchor", "middle")
                .style("fill","black")
                .attr("alignment-baseline","auto");   
            
        }
        if (cax.y){
            this.y_axis_svg=this.graph_area.append("g");
            this.y_scale = scaleLinear();
            this.y_axis_call= axisLeft(this.y_scale);
            this.y_axis_label=this.svg.append("text")
                .attr("text-anchor", "middle")
                .style("fill","black")  // 
                .attr("alignment-baseline","auto");  
     

        }
          
      
	}

    
    _setMargins(){
        const ax = this.config.axis
        this.margins = {
            top: ax.tx?ax.tx.size:10,
            right: ax.ry?ax.ry.size:10,
            bottom: ax.x?ax.x.size:10,
            left: ax.y?ax.y.size:10
        }
    }

    _getContentDimensions(){
           
        return {
            top:this.margins.top,
            left:this.margins.left,
            height:this.height-this.margins.bottom- this.margins.top,
            width:this.width-this.margins.left-this.margins.right
            
        }
    }

    updateAxis(){     
        const  dim = this._getContentDimensions();
        const ax = this.config.axis;
        if (ax.x){
            this.x_scale.range([0, dim.width]);
            this.x_axis_svg.attr("transform", `translate(0,${dim.height})`);
            this.x_axis_svg.transition().call(this.x_axis_call);
            this.x_axis_label.attr("transform",`translate(${this.width/2},${this.height-4})`)
                .text(this.config.axis.x.label)
                .attr("font-size",this.config.axis.x.textsize+"px")

        }

        if (ax.y){
            this.y_scale.range([0, dim.height]);
            this.y_axis_svg.transition().call(this.y_axis_call);
            this.y_axis_label.attr("transform", `translate(10,${dim.height/2}) rotate(-90)`)
                .text(this.config.axis.y.label)
                .attr("font-size",this.config.axis.y.textsize+"px")
    
        }
    }

    setSize(x,y){
        super.setSize(x,y);
        if (this.svg){
            this.svg 
            .attr("width", this.width)
            .attr("height", this.height)
        }
   
    }

    getImage(callback,type){
        if (this.addToImage){
            this.addToImage();
        }
        const cl = this.config.color_legend;
        //temporarily attach legend
        let leg_g =null
        let leg=null;
        if (this.legend){
            leg = this.legend.querySelector("g");
            leg_g =  select(leg).attr("transform",`translate(${cl.pos[0]},${cl.pos[1]})`);
            this.svg.node().append(leg_g.node());         
        }
		const svgAsXML = (new XMLSerializer).serializeToString(this.svg.node());

        if (type==="png"){
            this.getImageFromSVG(this.svg.node(),(canvas)=>{        
                callback(canvas);
                if (this.removeFromImage){
                    this.removeFromImage();
                }
                if (leg_g){
                    leg_g.attr("transform",`translate(0,0)`);
                    this.legend.querySelector("svg").append(leg);
                }
            })
        }
        else{
            callback("data:image/svg+xml," + encodeURIComponent(svgAsXML));
            if (this.removeFromImage){
                this.removeFromImage();
            }  
            if (leg_g){
                leg_g.attr("transform",`translate(0,0)`);
                this.legend.querySelector("svg").append(leg);
            }
        }    
    }

}

export {SVGChart};