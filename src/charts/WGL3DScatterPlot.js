import { WGL2DI } from "../webgl/WGL2DI.js";
import WGLChart from "./WGLChart.js";
import BaseChart from "./BaseChart.js";

class WGL3DScatterPlot extends WGLChart{
    constructor(dataStore,div,config){
        if (!config.title){
            config.title= dataStore.getColumnName(config.param[0])+" x "+
            dataStore.getColumnName(config.param[1])+" x "+
            dataStore.getColumnName(config.param[2]);
        }
		super(dataStore,div,config,{});
        const c = this.config;
         //get the x,y,z columns
        const p = c.param;
        this.dim = this.dataStore.getDimension("range_dimension",p);
        const ranges = p.map(x=> this.dataStore.getColumnRange(x))
        const max= Math.max(...ranges);
        this.defaultCDistance = max*4;
        c.brush = c.brush || "default";
		this.app= new WGL2DI(this.graphDiv,{mode:"3d",brush:c.brush,cameraDistance:this.defaultCDistance});
       

        this.app.addHandler("brush_stopped",(indexes,is_poly)=>{
		    this.resetButton.style.display = "inline";
			this.app.setFilter(true);
		    if (!is_poly){
		        this._createFilter(indexes);
		    }
		    else{
		    	this._createPolyFilter(range);
		    }
		});
       
        const colorFunc = this.afterAppCreation();

        this.app.addCircles({
            x:dataStore.getRawColumn(p[0]),
            y:dataStore.getRawColumn(p[1]),
            z:dataStore.getRawColumn(p[2]),
            localFilter:this.dim.getLocalFilter(),
            globalFilter:dataStore.getFilter(),
            colorFunc:colorFunc
          
          });

          c.radius=c.radius || 5;
          c.opacity= c.opacity || 0.8;
          
          this.app.setPointRadius(c.radius);
          this.app.setPointOpacity(c.opacity);
          c.axis_scales=c.axis_scales || [1,1,1];
          this.app.axisScales=c.axis_scales;

          this.addAxis();
          this.centerGraph(c.camera);
          this.app.refresh();
          
	}
    setAxisScales(index,val){
        this.config.axis_scales[index]=val;
        this.app.refresh();
    }

    getSettings(){        
        let settings =  super.getSettings({pointMax:30,pointMin:1})
        const c = this.config;
        settings = settings.concat([
            {
                type:"slider",
                max:3,
                min:0.5,
                doc:this.__doc__,
                current_value:c.axis_scales[0],
                label:"X Axis Scale",
                func:(x)=>{
                    this.setAxisScales(0,x)
                }
              
            },
            {
                type:"slider",
                max:3,
                min:0.5,
                doc:this.__doc__,
                current_value:c.axis_scales[1],
                label:"Y Axis Scale",
                func:(x)=>{
                    this.setAxisScales(1,x)
                }
            },
            {
                type:"slider",
                max:3,
                min:0.5,
                doc:this.__doc__,
                current_value:c.axis_scales[2],
                label:"Z Axis Scale",
                func:(x)=>{
                    this.setAxisScales(2,x)
                }
            }
        ]);
        return settings;   
    }

    centerGraph(values){
        if (values){
            this.app.setCamera(values.distance,values.theta,values.phi);
        }
        else{
            this.app.setCamera(this.defaultCDistance,-1.038,0.261);
        }
        
    }

    getConfig(){
        const config = super.getConfig();
        config.camera=this.app.getCameraSettings();
        return config;
    }

    _createFilter(indexes){
        this.dim.filter("filterOnIndex",[],indexes);
    }

    onDataAdded(newSize){
        const p = this.config.param;
		const config = this.getSetupConfig();
        config.x=this.dataStore.getRawColumn(p[0]);
        config.y=this.dataStore.getRawColumn(p[1]);
        config.z=this.dataStore.getRawColumn(p[2]);
        this.app.updateSize(newSize,config);
        super.onDataAdded(newSize);
    }

    addAxis(){
        for (let index=0;index<3;index++){
            let mm = this.dataStore.getMinMaxForColumn(this.config.param[index]);
            const from = [0,0,0];
            from[index]=mm[0];
            const to = [0,0,0];
            to[index]=mm[1];
            const color = [0,0,0];
            color[index]=255;
            this.app.addLine(from,to,color);
        }   
    }
}

BaseChart.types["wgl_3d_scatter_plot"]={
    name:"3D Scatter Plot",
    class:WGL3DScatterPlot,
    params:[{
        type:"number",
        name:"X axis"
    },
    {
        type:"number",
        name:"Y axis"
    },
    {
        type:"number",
        name:"Z axis"
    },

    ]
}

export default WGL3DScatterPlot;