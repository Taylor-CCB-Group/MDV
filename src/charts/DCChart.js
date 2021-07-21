import {BaseChart} from "./BaseChart.js";
import * as dc from "dc";

class DCChart extends BaseChart{
	constructor(dataStore,chartType,div,config){
    	super(dataStore,div,config);
        let self=this;
        this.chart=chartType(this.contentDiv);
        this.chart.controlsUseVisibility(true);
        this.chart.on("filtered",()=>{
            this.resetButton.style.visibility=this.chart.hasFilter()?"visible":"hidden"
        })
    }

    clearFilters(){
        this.chart.filterAll();
        dc.redrawAll();
    }
    


    setSize(x,y){
        super.setSize(x,y);
        this.chart.width(this.width).height(this.height);    
    }


     /** Removes the chart from the dom, clearing any filters
    * The updated items are returned
    * 
    */
    remove(){   
        this.dim.filter(null);
        let left = this.dim.getIds();
        this.dim.dispose();
        this.chart.resetSvg();
        this.div.remove();
        return left;
    }

    dataChanged(not_broadcast){
        let self = this;
        this.not_broadcast=true;

        let filter= this.chart.filters();
          if (filter.length>0){
           this.chart.filter(null);
        }
        
		this.setParameters();
		this.not_broadcast=false;
      
        if (filter.length>0){
             this.not_broadcast=not_broadcast
           
            this.chart.filter([filter]);
        }
    }


    getCanvas(callback){
    	let sc =null;
    	let svg = this.chart.svg().node();
		if (this.config.color_by){
		    sc = "#"+this.getColorOptions().div+"-bar";
			svg.append(d3.select(sc).attr("x",70).node());
			
		}
    	
        let self = this;
    	getImageFromSVG(svg,function(canvas,ctx){
    	   callback(canvas);
    	    if (sc){
        	$("#"+self.getColorOptions().div).append($(sc));
        }

    	});

    }

    getSVG(callback){  	
    	let sc =null;
    	let svg = this.chart.svg().node();
		if (this.config.color_by){
		    sc = "#"+this.getColorOptions().div+"-bar";			
			svg.append(d3.select(sc).attr("x",70).node());
			
		}
	
		 let svgAsXML = (new XMLSerializer).serializeToString(svg); 
         callback("data:image/svg+xml," + encodeURIComponent(svgAsXML));
         if (sc){
        	$("#"+this.getColorOptions().div).append($(sc));
        }
    }
    
	colorByField(params){
		if (!params){
			delete this.config.color_by;
			let div = this.getColorOptions().div;
			//$("#"+div+"-bar").remove();
		}
		else{
			this.config.color_by={
    			column:params.column,
    			scheme:params.scheme,
    			value_to_color:params.value_to_color
    		}
		}

    	this.setParameters();
	}
}

export {DCChart};
