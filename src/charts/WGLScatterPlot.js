
import {WGL2DI} from "../webgl/WGL2DI.js";
import {WGLChart} from "./WGLChart.js";
import {BaseChart} from "./BaseChart.js"


class WGLScatterPlot extends WGLChart{
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
	

    
        this.x= this.config.param[0];
        this.y= this.config.param[1];
        this.dim = this.dataStore.getDimension("range_dimension",this.config.param); 
        this.minMaxX=this.dataStore.getMinMaxForColumn(this.x);
        this.minMaxY=this.dataStore.getMinMaxForColumn(this.y);
		this.type="wgl_scatter_plot";

		let self = this;
		const c = this.config;
		c.brush = c.brush || "poly";

        let appConf= {brush:c.brush};
        
        this.app = new WGL2DI(this.graphDiv,appConf);
            this.app.addHandler("zoom_stopped",function(data){
            self._updateScale(data);
			self.updateAxis();
          
        });

		const colorFunc = this.afterAppCreation();

        this.app.addHandler("panning_stopped",function(data){
            self._updateScale(data);
			self.updateAxis();
          
        });
	
        this.app.addHandler("brush_stopped",function(range,is_poly){
		    self.resetButton.style.display = "inline";
			self.app.setFilter(true);
		    if (!is_poly){
		        self._createFilter(range);
		    }
		    else{
		    	self._createPolyFilter(range);
		    }
		});


		if (this.config.background_image){
			this.addBackgroundImage();
		}

		
	

		this.app.addCircles({
			x:this.dataStore.getRawColumn(this.x),
			y:this.dataStore.getRawColumn(this.y),
			localFilter:this.dim.getLocalFilter(),
			globalFilter:this.dataStore.getFilter(),
			colorFunc:colorFunc
			
		});
		this.defaultRadius=this._calculateRadius();
		console.log(this.defaultRadius)
	
		this.config.radius= 5;//this.config.radius || this.defaultRadius;
		this.config.opacity= this.config.opacity || 0.8;
		
		this.app.setPointRadius(this.config.radius);
		this.app.setPointOpacity(this.config.opacity);
		
		this.centerGraph();
		this.app.refresh();
		   
    }




	addBackgroundImage(){
		const ic = this.config.background_image;
		const im = new Image();
		im.src= ic.url;
		im.onload=()=>{
			this.app.addImage(im,ic);
			this.app.refresh();
		}
	}




    highlightDataItem(key){
    	this.app.highlightPoint(key);
    }

    getFilter(){
    	let filters=[];
    	if (!this.range){
    		return filters;
    	}
    	filters.push({field:this.config.param[0],operand:"between",value:[this.range.x_min,this.range.x_max]});
    	filters.push({field:this.config.param[1],operand:"between",value:[this.range.y_min,this.range.y_max]});
    	return filters;
    }





    _updateScale(range){	
        this.x_scale.domain([range.x_range[0],range.x_range[1]]);     
        this.y_scale.domain([-range.y_range[0],-range.y_range[1]]);   
    }



    _calculateRadius(){
    	let width = this.width?this.width:1;
		let max_x=this.config.max_x || this.config.max_x==0?this.config.max_x:this.minMaxX[1]; 	
    	let min_x=this.config.min_x || this.config.min_x==0?this.config.min_x:this.minMaxX[0];
    	if (this.config.axis.x_log_scale){
    		max_x= this._getLogValue(max_x);
    		min_x= this._getLogValue(min_x);
    	}
		
		let range = max_x-min_x
		let pt_px = Math.pow(range,1.50)/(Math.pow(this.dataStore.size,0.8)+range);
		//pt_px = pt_px<0.001?0.001:pt_px;
		return (pt_px)*5;
		
    }

  

    removeFilter(){
    	if (this.filter===null){
    		return;
    	}
		this.dim.removeFilter()	
    	this.app.clearBrush();
        this.app.setFilter(false);
		this.app.refresh();
    }
    
    _createFilter(range){
		this.range=range;
	
	
        if (range==null){
            this.dim.removeFilter();
        }
        else{
        	let y_max=range.y_max;
        	let y_min=range.y_min;
        	let x_max=range.x_max;
        	let x_min=range.x_min;


         	if (this.config.axis.y_log_scale){
         		y_max=this._getInverseLogValue(y_max);
         		y_min=this._getInverseLogValue(y_min);
         	}
         	if (this.config.axis.x_log_scale){
         		x_max=this._getInverseLogValue(x_max);
         		x_min=this._getInverseLogValue(x_min);
         	}
         	
        
	
			this.filter= [[x_min,x_max],[y_min,y_max]]
            this.dim.filterSquare(this.filter[0],this.filter[1]);
		   

        }
    }



    _createPolyFilter(vs){
    	this.range=true;
    	for (let pt of vs){
    		pt[1]=-pt[1];
    		if (this.config.axis.x_log_scale){
    			pt[0]=this._getInverseLogValue(pt[0]);
    		}
    		if (this.config.axis.y_log_scale){
    			pt[1]=this._getInverseLogValue(pt[1]);
    		}
    	}
		this.filter= vs;
		this.dim.filterPoly(vs);
    }



    centerGraph(){
    	let max_x=this.config.max_x || this.config.max_x==0?this.config.max_x:this.minMaxX[1];
    	let max_y=this.config.max_y || this.config.max_y==0?this.config.max_y:this.minMaxY[1];
    	let min_x=this.config.min_x || this.config.min_x==0?this.config.min_x:this.minMaxX[0];
    	let min_y=this.config.min_y || this.config.min_y==0?this.config.min_y:this.minMaxY[0];
    	if (this.config.axis.y_log_scale){
    		max_y=this._getLogValue(max_y);   			   		
    		min_y=this._getLogValue(min_y);    		
    	}
    	if (this.config.axis.x_log_scale){
    		max_x=this._getLogValue(max_x);   		
    		min_x=this._getLogValue(min_x);
    	}

        let x_margin=((max_x-min_x)/20);
        let y_margin=((max_y-min_y)/20);
        let x_range = (max_x-min_x)+2*x_margin;
        let y_range= (max_y-min_y)+2*y_margin;

        const dim = this._getContentDimensions();

        this.app.x_scale =(dim.width)/x_range;
        this.app.y_scale =(dim.height)/y_range;
        this.app.offset[0]=-(min_x-x_margin);
        this.app.offset[1]=(max_y+y_margin);
		this._updateScale(this.app.getRange());
		this.updateAxis();

    }

	getSettings(){
        const settings = super.getSettings({pointMax:30,pointMin:0});
        const c = this.config;
		return settings.concat([
			{
				type:"radiobuttons",
				label:"Brush Type",
				choices:[["Free Draw","poly"],["Rectangle","default"]],
				current_value:c.brush,
				func:(x)=>{
					this.app.clearBrush();
					this.app.config.brush=x;
					c.brush=x;
				}

			}
		]);
	}

	
}

BaseChart.types["wgl_scatter_plot"]={
    name:"2D Scatter Plot",
    class:WGLScatterPlot,
    params:[{
        type:"number",
        name:"X axis"
    },
    {
        type:"number",
        name:"Y axis"
    }
    ]
}


export {WGLScatterPlot}