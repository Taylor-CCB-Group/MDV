import {SVGChart} from "./SVGChart.js";
import {createEl} from "../utilities/Elements.js";
import { hexToRGB } from "../datastore/DataStore.js";
import {select} from "d3-selection";
import "d3-transition";
class WGLChart extends SVGChart{
    constructor(dataStore,div,config,axisTypes){
        super(dataStore,div,config,axisTypes);
        const c  = this.config;
		if (!c.tooltip){
			c.tooltip={show:false}
		}

		c.default_color= c.default_color || "#377eb8";
          
        
        const box =this._getContentDimensions();
		this.graphDiv = createEl("div",{
            styles:{
                position:"absolute",
				left:box.left+"px",
				top:box.top+"px",
				width:box.width+"px",
				height:box.height+"px"
            }
        },this.contentDiv);
        this.tooltip = createEl("div",{
			styles:{
				position:"absolute",
                zIndex:100,
				display:"none",
				"background-color":"lightgray"
			}
		},this.graphDiv)

    }
    
    afterAppCreation(){
        const self = this;
        const c= this.config;
        this.app.addHandler("object_over",function(e,index){
			if (c.tooltip.show){
				self.showTooltip(e,index)
			}       
        });
		this.app.addHandler("object_out",function(e,index){
            self.tooltip.style.display="none"        
        });

        c.default_color= c.default_color || "#377eb8";

        const dColor = hexToRGB(c.default_color)
		let colorFunc=()=>dColor;
		
		if (c.color_by){
			colorFunc=this.dataStore.getColorFunction(c.color_by,{asArray:true});
            if (!c.color_legend){
                c.color_legend={
                    display:true
                }
            }
            if (c.color_legend.display){
                this.addColorLegend(c.color_by);
            }
		}

      
        return colorFunc;

    }

    showTooltip(e,index){
		var rect = this.graphDiv.getBoundingClientRect();
		const row =  this.dataStore.getRowAsObject(index);
    	const x = e.clientX-rect.left;
		const y =e.clientY-rect.top;
		this.tooltip.style.display="block";
		this.tooltip.style.left = x+"px";
		this.tooltip.style.top = y+"px";
		this.tooltip.textContent=row[this.config.tooltip.column];

	}

    setSize(x,y){
		super.setSize(x,y);
		const dim = this._getContentDimensions();
		this.app.setSize(dim.width,dim.height);
		this.updateAxis();
	}

    onDataFiltered(dim){
        if (dim === "all_removed"){
            this.app.clearBrush();
            this.app.setFilter(false);
            this.resetButton.style.display="none";
            
        }      
        this.app.refresh();		
	}

    remove(){    
        this.dim.destroy();
        this.app.remove();
		super.remove()       
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


    addColorLegend(){
        if (!this.config.color_legend.display){
            if (this.legend){
                this.legend.remove();
            }
            return;
        }
        const box = this._getContentDimensions();
        let lt= 0;
        let ll = 0;
        if (this.legend){
            ll = this.legend.style.left;
            lt = this.legend.style.top;
            this.legend.remove();
        }
        else {
            const cl= this.config.color_legend;
            if (!cl.pos){
                cl.pos=[box.left,box.top]
            }
            ll= cl.pos[0]+"px";
            lt= cl.pos[1]+"px";
        }
        this.legend = this.dataStore.getColorLegend(this.config.color_by);
        this.contentDiv.append(this.legend);
       
        this.legend.style.left= ll;
        this.legend.style.top= lt;
        this.legend.__doc__=this.__doc__;

    }

    colorByColumn(column){
		const colorFunc = this.dataStore.getColorFunction(column,{asArray:true});
        if (!this.config.color_legend){
            this.config.color_legend={
                display:true
            }
        }
        
        this.addColorLegend();
        
        
		const t = performance.now();
		this.app.colorPoints(colorFunc);
		console.log("color time:"+(performance.now()-t))
        setTimeout(()=>{
            this.app.refresh();
        },50);
        
	}

	colorByDefault(){
        if (this.legend){
            this.legend.remove();
        }
		const col = hexToRGB(this.config.default_color);
		this.app.colorPoints(()=>col);
		setTimeout(()=>{
            this.app.refresh();
        },50);

	}

    changeBaseDocument(doc){
		this.app.clearBrush();
        this.app.__doc__=doc;
		super.changeBaseDocument(doc);
	}
    
    getColorOptions(){		
		return {
			colorby:"all",
			has_default_color:true
		}
	}
    
    pinChart(){
        this.isPinned=true;
        //store a copy of the global filter and pass it
        //to the shader
        this.pinnedFilter=this.dataStore.getFilter().slice(0);
        this.app.setGlobalFilter(this.pinnedFilter);
    }

    unpinChart(){
        this.isPinned =false;
        this.pinnedFilter=null;
        this.app.setGlobalFilter(this.dataStore.getFilter());
        this.app.refresh();
    }

    addToImage(){
        this.app.refresh();
    	let url =this.app.canvas.toDataURL('image/png');
		url = url.replace(/^data:image\/png/,'data:application/octet-stream');
		const box = this._getContentDimensions();
        
		this.svg.append('image')
		    .attr("xlink:href",url)
		    .attr("width",box.width)
		    .attr("height",box.height)
		    .attr("x",box.left)
		    .attr("y",box.top);
    }
    removeFromImage(){
        this.svg.select("image").remove();
    }

   
    getSettings(conf){
        const settings = super.getSettings();
        const c = this.config;
        const cols = this.dataStore.getColumnList();
    
        return settings.concat([
            {
                type:"slider",
                max:1,
                min:0,
                doc:this.__doc__,
                current_value:c.opacity,
                label:"Point Opacity",
                func:(x)=>{
                    c.opacity=x;
                    this.app.setPointOpacity(x)
                    this.app.refresh();
                }
            },
            {
                type:"slider",
                max:conf.pointMax,
                min:conf.pointMin,
        
                doc:this.__doc__,
                current_value:c.radius,
                label:"Point Size",
                func:(x)=>{
                    c.radius=x;
                    this.app.setPointRadius(x)
                    this.app.refresh();
                    console.log(x);
                }
            },
            {
                type:"check",
                label:"Show Tooltip",
                current_value:c.tooltip.show,
                func:(x)=>{
                    c.tooltip.show=x;
                    if (!c.tooltip.column){
                        c.tooltip.column= cols[0].field;
                    }
                }
            },
            {
                type:"dropdown",
                label:"Tooltip value",
                current_value:c.tooltip.column || cols[0].field,
                values:[cols,"name","field"],
                func:(x)=>{
                    c.tooltip.column=x;
                }

            },
          
            {
                type:"button",
                label:"Centre Plot",
                func:(x)=>{
                    this.centerGraph();
                    this.app.refresh();
                }

            }


        ]);
    }
    
}




export {WGLChart}