import {ColorStackChart} from "./ColorStackChart.js";
import * as dc from "dc";
import * as d3 from "d3";


class HistogramChart1 extends ColorStackChart{
    constructor(dataStore,div,config){    	 
        super(dataStore,dc.barChart,div,config);     
    }

    init(){
    	this.type="bar_chart";
        let self =this;
        
        if (this.config.bin_number === undefined){
        	this.config.bin_number=10;
        }
        

		let min_max= this.dataStore.getMinMaxForColumn(this.config.param);
		this.min = min_max[0];
		if (this.config.display_min === undefined){
			this.config.display_min=min_max[0];
		}
		this.max = min_max[1];
		if (this.config.display_max === undefined){
			this.config.display_max = min_max[1];
		}

    }

   
    setParameters(params){
        let self=this;
        if (params){
        	this._setParam(params,"display_max");
         	this._setParam(params,"display_min");       	
			this._setParam(params,"bin_number");
			this._setParam(params,"max_y");
			
     
        }
        if (this.dim){
        	this.dim.dispose();
        }
         
        this.range = this.config.display_max-this.config.display_min

        this.bin_width=(this.range/(this.config.bin_number));
        let fudge = this.bin_width/10;
        this.dim=this.dataStore.getSimpleDimension(this.config.param,[this.config.display_min,this.config.display_max])
      
        if (this.group){
            this.group.dispose()
        }
        this.group = this.dim.group(function(d){
            let v = (self.bin_width * Math.floor(d/self.bin_width));
            return v;

        });
        
        this.chart.dimension(this.dim)
                   .xUnits(dc.units.fp.precision(this.bin_width))
                   .x(d3.scaleLinear().domain([this.config.display_min-this.bin_width,this.config.display_max+this.bin_width]))
                   
                   .yAxisLabel("",0)
                    .xAxis().ticks(Math.round(this.width/30)).tickFormat(function(v,i){
                        if (Math.abs(v)>=1000){
                        if ((i+1)%2==0){
                            if (Math.abs(v)>=100000){
                                return Number.parseFloat(v).toPrecision(2);
                            }
                            return v;
                        }
                        return "";
                        }
                        else{
                            return v
                        }
                    });

        
        this.chart.yAxis().ticks(Math.round(this.height/25)).tickFormat(function(v,i){
                            if (v>=100000){
                                return Number.parseFloat(v).toPrecision(2);
                            }
                            return v;
                        
         
                    });

      
        this.chart.margins().right = 10;
        this.chart.margins().left=40;
        this.setColorStack();
        this.chart.render();     
         
    }
    setSize(x,y){
        super.setSize(x,y);
         this.chart.x(d3.scaleLinear().domain([this.config.display_min-this.bin_width,this.config.display_max+this.bin_width]))
         .xAxis().ticks(Math.round(this.width/30));
         this.chart.yAxis().ticks(Math.round(this.height/25));
         this.chart.redraw();
    }
}

export {HistogramChart1}
