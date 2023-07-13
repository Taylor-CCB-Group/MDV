

import { createEl } from "../utilities/Elements";
import BaseChart from "./BaseChart";
import SVGChart from "./SVGChart";
import { WGL2DI } from "../webgl/WGL2DI";


const color_scheme=["#313695","#4575B4","#74ADD1","#ABD9E9", "#E0F3F8", "#E0F3F8","#FFFFBF", "#FEE090","#FDAE61","#F46D43" ,"#D73027" ,"#A50026" ].reverse();

class DeepToolsHeatMap extends SVGChart{
    constructor(dataStore,div,config){
        super(dataStore,div,config,{x:{type:"band"},y:{custom:true}});
        const box =this._getContentDimensions();
        //create the div to house the webgl
        this.graphDiv = createEl("div",{
            styles:{
                position:"absolute",
                left:box.left+"px",
                top:box.top+"px",
                width:box.width+"px",
                height:box.height+"px"
            }
        },this.contentDiv);
        const c = this.config;
        if (!c.color_scale){
            c.color_scale={
                range:[0,255],
                log:false            
            }
        }
        //create the webgl
        this.app = new WGL2DI(this.graphDiv,{lock_x_axis:true});

        this.app.addHandler("object_clicked",i=>{
            //get the position in flat non sparse array array
            const p = this.positions[i];
            // work out the row
            const r = Math.floor(p/this.hm_config.cols);
            //get the locus index
            const index = this.map_index[r];
            this.dataStore.dataHighlighted([index],this);
            
          //wgl.offset=[10,-r*20];
          //wgl.y_scale=wgl.height/(10*20);
          //wgl.addLine([0,r*20],[cols*20+600,r*20],[255,255,0]);
          //wgl.addLine([0,r*20+20],[cols*20+600,r*20+20],[255,255,0])
          //wgl.refresh();
          
   });
        //load the correct data
        this.hm_config =this.dataStore.deeptools.maps[c.heatmap];
        this.x_scale.domain(this.hm_config.groups.slice(0));
       
        this.dataStore.loadBinaryData(this.hm_config.data).then(buff=>{
            const r = this.hm_config.rows;
            const c =  this.hm_config.cols;
            this.map_index = new Uint32Array(buff,0,r);
            this.reverse_map= new Map();
            for (let n=0;n<this.map_index.length;n++){
                this.reverse_map.set(this.map_index[n],n)
            }
            /*
            const dd = this.dataStore.getRawColumn("TSS_RE");
            const vs = this.dataStore.getColumnValues("TSS_RE");
            const mi = this.map_index;
            let gr= vs[dd[mi[0]]]
            let st=0
            const grps=[{"group":gr,"start":0}]
            for (let i=0;i<mi.length; i++){
                const a=dd[mi[i]];
                if (mi[i]!==0 && vs[a] !==gr){
                    console.log(i+":"+gr+":"+vs[a]);
                    grps.push({"group":vs[a],"start":i});
                    gr=vs[a];
                 
                }
            }
            */
            this.data= new Uint8Array(buff,r*4,r*c);
            const len = this.data.length;
            let z=0;
            for (let n=0;n<this.data.length;n++){
                if (this.data[n]===0){
                    z++
                }
            }
            this.colors= new Uint8Array(z*3)
            this.x= new Float32Array(z);
            this.y= new Float32Array(z);
            this.positions= new Float32Array(z);
            const cpg = c/this.hm_config.groups.length;
            let i=0;
            for (let n=0;n<len;n++){
                if (this.data[n]===0){
                    continue;
                }
                const col = n%c;
                const row = Math.floor(n/c);
                this.x[i]=col*20 + (Math.floor(col/cpg)*200)+10;
                this.y[i]=row*20+10;
                this.positions[i]=n;
                i++;
            }
            this.setColor();
            const dim = this._getContentDimensions();
            this.app.x_scale= dim.width/((c*20)+(this.hm_config.groups.length-1)*200);
            this.app.y_scale = dim.height/(r*20);
            const gl = this.hm_config.groups.length;
            const gw= (c/gl)*20;
            let x=0;
            const base_color = this.colorScale(0); 
            for (let n=0;n<gl;n++){
                this.app.addRectangle([x,0],gw,r*20,base_color);
                x+=gw+200;
            }
            
           
           
            this.app.addSquares({
                x:this.x,
                y:this.y,
                colors:this.colors
            });
            this.app.refresh()              
        });
    }

    setSize(x,y){
		super.setSize(x,y);
        //rescale the x axis
        //const rs = (w,h)=>[w/((this.hm_config.cols*20)+(this.hm_config.groups.length-1)*200),h/(this.hm_config.rows*20)]
		const dim = this._getContentDimensions();
		this.app.setSize(dim.width,dim.height);
        this.graphDiv.style.left = dim.left+"px";
        this.graphDiv.style.top = dim.top+"px";
		this.updateAxis();
	}


    onDataHighlighted(data){
        this.app.removeAllLines();
        const r = this.reverse_map.get(data.indexes[0]);
        const cols= this.hm_config.cols;
        if (r){
            this.app.addLine([0,r*20],[cols*20+600,r*20],[255,255,0]);
            this.app.addLine([0,r*20+20],[cols*20+600,r*20+20],[255,255,0])
        }
        
        if (data.source!==this && r){
            this.app.offset=[0,-r*20 +300];
            this.app.y_scale=this.app.height/(30*20);
        }
        this.app.refresh();

    }


    setColor(){
        const r= this.config.color_scale.range;
        const conf ={
            datatype:"integer",
            asArray:true,
            overideValues:{
                colorLogScale:this.config.color_scale.log,
                colors:color_scheme,
                min:r[0],
                max:r[1],        
                bins:200,      
            }
        }
        this.colorScale = this.dataStore.getColorFunction(null,conf);
        const c = this.colors;
        for (let n=0;n<this.data.length;n++){
            const i =this.positions[n];
            const col = this.colorScale(this.data[i]);
            const p = 3*n
            c[p]=col[0];
            c[p+1]=col[1];
            c[p+2]=col[2]
        }
    }


    getSettings(){
        let settings= super.getSettings();
        const c = this.config
        const cs= c.color_scale;
        settings = settings.concat([{        
            type:"doubleslider",
            max:255,
            min:0,
            doc:this.__doc__,
            current_value:cs.range,
            label:"Color Scale",
            func:(x,y)=>{     
                cs.range=[x,y];
                this.setColor();
                this.app.refresh();
            }
        }]);
        return settings;
    }
}

BaseChart.types["deeptools_heatmap"]={
    class:DeepToolsHeatMap,
    required:["deeptools"],
    name:"DeepTools HeatMap",
    extra_controls:(dataSource)=>{
        const values = Object.keys(dataSource.deeptools.maps).map(x=>({name:x,value:x}))
        return [
            {
                type:"dropdown",
                name:"map",
                label:"Map",
                values
            }
        ];
    },

    init:(config,dataSource,ec)=>{
        config.heatmap= ec["map"]
    }
}
