
import {WGL2DI} from "../webgl/WGL2DI.js";
import WGLChart from "./WGLChart.js";
import BaseChart from "./BaseChart.js"


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
        this.dim = this.getDimension()
		const bf = config.background_filter;
		if (bf){
			this.dim.setBackgroundFilter(bf.column,bf.category);
		}
        this.minMaxX=this.dataStore.getMinMaxForColumn(this.x);
        this.minMaxY=this.dataStore.getMinMaxForColumn(this.y);
		this.type="wgl_scatter_plot";

		let self = this;
		const c = this.config;
		c.brush = c.brush || "poly";

        let appConf= {brush:c.brush};
        
        this.app = new WGL2DI(this.graphDiv,appConf);

        this.app.addHandler("zoom_stopped",data=>this.handlePanZoom(data));
		this.app.addHandler("panning_stopped",data=>this.handlePanZoom(data));
	
		this.centerGraph();
          
		const colorFunc = this.afterAppCreation();


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

		if (c.background_image){
			c.image_opacity = c.image_opacity==null?1:c.image_opacity;
			this.addBackgroundImage(c.background_image);
		}

		this.app.addCircles({
			x:this.dataStore.getRawColumn(this.x),
			y:this.dataStore.getRawColumn(this.y),
			localFilter:this.dim.getLocalFilter(),
			globalFilter:this.dataStore.getFilter(),
			colorFunc:colorFunc
			
		});

		if (c.offsets){
			const col = this.dataStore.columnIndex[c.offsets.param];
			const v = col.values.indexOf(c.offsets.value)
			this.app.setOffsets(col.data,v,c.offsets.offsets);
		}
		this.defaultRadius=this._calculateRadius();
		console.log(this.defaultRadius);
	
		c.radius= c.radius || 2;
		
		c.opacity= c.opacity==null?0.8:c.opacity
		
		this.app.setPointRadius(this.config.radius);
		this.app.setPointOpacity(this.config.opacity);
		
		this.onDataFiltered();	   
    }

	getDimension(){
		return this.dataStore.getDimension("range_dimension");
	}

	handlePanZoom(range){
		if (range.imageMoved){
			const i = range.imageMoved;
			const c=this.config
			const name =c.background_image.name;
			const inf = c.image_choices.find(x=>x[0]===name);
			i.position[1]=-(i.position[1]+i.height)
			c.background_image.height=i.height;
			c.background_image.width= i.width;
			c.background_image.position= i.position;

			inf[2]=i.width;
			inf[3]=i.height;
			inf[4]=i.position[0];
			inf[5]=i.position[1];

		}
		else{
			this._updateScale(range);
			this.updateAxis();
		}
	}

	onDataAdded(newSize){
		const config = this.getSetupConfig();
		config.x=this.dataStore.getRawColumn(this.x);
		config.y=this.dataStore.getRawColumn(this.y);
        this.app.updateSize(newSize,config);
        super.onDataAdded(newSize);
    }

	addBackgroundImage(ic,change=false){
		
		const im = new Image();
		const c = this.config;
		im.src= ic.url;
		im.onload=()=>{
			if (change){
				this.app.changeImage(im,ic,0);
			}
			else{
				this.app.addImage(im,ic);
				this.app.changeImageOpacity(0,c.image_opacity);
			}		
			this.app.refresh();
		}
	}

	linkToOtherChart(chart){
		this.app.addHandler("pan_or_zoom",(offset,x_scale,y_scale)=>{
			chart.app.offset=[offset[0],offset[1]];
			chart.app.x_scale=x_scale;
			chart.app.y_scale=y_scale;
			chart.app.refresh();
			chart._updateScale(this.app.getRange());
			chart.updateAxis();
		},this.config.id)
	}



    highlightDataItem(key){
    	this.app.highlightPoint(key);
    }

    getFilter(){
    	
    	if (!this.range || this.range === true){
    		return null
    	}
		const f= {}
    	f[this.config.param[0]]=[this.range.x_min,this.range.x_max];
    	f[this.config.param[1]]=[this.range.y_min,this.range.y_max];
    	return f;
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
            this.dim.filter("filterSquare",this.config.param,{range1:this.filter[0],range2:this.filter[1]});
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
		this.dim.filter("filterPoly",this.config.param,vs);
    }



    centerGraph(){
		const roi = this.config.roi || {};
    	let max_x=roi.max_x || this.minMaxX[1];
    	let max_y=roi.max_y || this.minMaxY[1];
    	let min_x=roi.min_x || this.minMaxX[0];
    	let min_y=roi.min_y || this.minMaxY[0];
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
		this.app.pointScale= (this.app.x_scale+this.app.y_scale)/2;
        this.app.offset[0]=-(min_x-x_margin);
        this.app.offset[1]=(max_y+y_margin);
		this._updateScale(this.app.getRange());
		this.updateAxis();

    }



	_addSliderToSettings(settings,im,index){
		settings.splice(1,0,{
			type:"slider",
			current_value:im.opacity,
			min:0,
			max:1,
			func:x=>{
				im.opacity=x;
				this.app.changeImageOpacity(index,x)
				this.app.refresh();

			}
		})

	}

	remove(){
		if (this.config.offsets){
			const x = this.dataStore.getRawColumn(this.x);
			const y= this.dataStore.getRawColumn(this.y);
			const lf = this.dim.getLocalFilter();
			const o = this.app.offsets.offsets;
			const arr  = this.app.offsets.arr;
			const v = this.app.offsets.value;
			if (o[0]!== 0 && o[1] !==0){
				for (let n=0;n<lf.length;n++){
					if (arr[n]===v && lf[n]!==2){
						x[n]-=o[0];
						y[n]-=o[1];
					}

				}
			}
		}
		super.remove();
	}

	getConfig(){
		const c = super.getConfig();
		if (c.offsets){
			c.offsets.offsets = this.app.offsets.offsets;
		}
		return c;
	}

	getSettings(){
        const settings = super.getSettings({pointMax:30,pointMin:0});
        const c = this.config;
		if (c.image_choices){
			const ic= c.image_choices.slice(0);
			ic.unshift(["None","__none__"]);
			let cv  = "__none__";
			if (c.background_image){
				cv = c.background_image.url;
			}
			settings.splice(1,0,{
				type:"dropdown",
				current_value:cv,
				values:[ic,0,1],
				label:"Change Image",
				func:x=>{
					if (x==="__none__"){
						delete c.background_image;
						this.app.removeImages();
						this.app.refresh();
						return;
					}
					let replace = true;
					if (!c.background_image){
						c.background_image={};
						replace=false;
					}
					
					c.background_image.url=x;
					let t  = c.title.split("-")[0]
					const f= c.image_choices.find(i=>i[1]===x);
					c.background_image.name=f[0];
					if (f[2]){
						c.background_image.width=f[2];
						c.background_image.height=f[3];
						c.background_image.position=[f[4],f[5]];
					}
					
					this.setTitle(t+"-"+f[0])
					this.addBackgroundImage(c.background_image,replace);

				}
				

				
			})
		}
		if (c.background_image){
			settings.splice(2,0,{
				label:"Image Opacity",
				      type:"slider",
				       current_value:c.image_opacity,
				       max:1,
				       min:0,
				       func:x=>{
				           c.image_opacity=x;
				           this.app.changeImageOpacity(0,x);
						   this.app.refresh();
				       } 
					      
			})

		}		
			
			
		
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


export default WGLScatterPlot;