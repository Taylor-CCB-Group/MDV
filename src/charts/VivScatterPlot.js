import WGLScatterPlot from "./WGLScatterPlot.js";
import DensityScatterPlot from "./DensityScatterPlot.js";
import BaseChart from "./BaseChart.js";
import {createEl} from "../utilities/Elements.js";
import {RGBToHex,hexToRGB} from "../datastore/DataStore.js";
import {BaseDialog} from "../utilities/Dialog.js";
import noUiSlider from "nouislider";

export class ColorChannelDialog extends BaseDialog{
    constructor(viv,doc){
        const config={
            width:500,
            title:"Channels",
            doc:doc || document
        }
        super(config,viv);
    }
    init(viv){
        this.viv=viv;
        this.mainDiv=createEl("div",{},this.dialog)
        const channels  =this.viv.getChannels();
        for (let c of channels){
            this.addSlider(c);
        }

        const addDiv = createEl("div",{
            styles:{
                padding:"4px"
            }
        },this.dialog);
        const sel = createEl("select",{},addDiv);
        const ac = viv.getAllChannels();
        for (let n=0;n<ac.length;n++){
            createEl("option",{
                text:ac[n].Name ?? ac[n].ID,
                value:n
            },sel)
        }
        const cca= createEl("input",{
            styles:{
                width:"25px",
                height:"20px",
                padding:"0px",
                margin:"0px 4px"
            },
            type:"color",
            value:"#ff0000"
        },addDiv);
        createEl("button",{
            classes:["ciview-button-sm"],
            text:"Add Channel"
        },addDiv).addEventListener("click",()=>{
            const ch = this.viv.addChannel({index:parseInt(sel.value),color:cca.value});
            this.addSlider(ch);
        })

        
    }


    addSlider(item){
        const cont = createEl("div",{
            styles:{
                display:"flex",
                padding:"4px",
                alignItems: "center"
            }
        },this.mainDiv);
        createEl("span",{
            text:item.name,
            classes:["mdv-flex-dynamic"],
            styles:{flexBasis:"30%"}
            
        },cont)
        const sl = createEl("div",{
            classes:["mdv-flex-dynamic"],
            styles:{
                padding:"0px 10px",
                flexBasis:"70%",
                overflow:"visible"
            }
        },cont);
        const d = item.domains;
        noUiSlider.create(sl, {
            start: [item.contrastLimits[0],item.contrastLimits[1]],
            range: {
                'min': d ? d[0] : 0,
                'max': d ? d[1] : 300
            },
            step:1,
            document:this.config.doc,
            tooltips:true
        },cont);
        sl.noUiSlider.on("update",(values)=>{
            item.contrastLimits=[parseFloat(values[0]),parseFloat(values[1])];
            this.viv.setChannel(item);
        });
        const cc= createEl("input",{
            classes:["mdv-flex-fixed"],
            styles:{
                width:"25px",
                height:"20px",
                padding:"0px",
            },
            type:"color",
            value:item.color
        },cont);

        cc.addEventListener("input",()=>{
            item.color = cc.value;
            this.viv.setChannel(item);
        });
        const dc= createEl("input",{
            classes:["mdv-flex-fixed","mdv-checkbox"],
            type:"checkbox"
        },cont);
        dc.checked= item.channelsVisible;
        dc.addEventListener("click",()=>{
            item.channelsVisible= dc.checked;
            this.viv.setChannel(item);
        })
        const del = createEl("i",{
            classes:["mdv-flex-fixed","fas","fa-times"],
            
            styles:{
                width:"16px",
                margin:"0px 3px",
                cursor:"pointer"
            }
        },cont);
       
        del.addEventListener("click",()=>{
            this.viv.removeChannel(item);
            cont.remove();

        })


    }
}


class VivScatterPlot extends DensityScatterPlot{
    constructor(dataStore,div,config){
        const x_name= dataStore.getColumnName(config.param[0]);
        const y_name = dataStore.getColumnName(config.param[1]);
        if (!config.axis){

            config.axis={
                x:{size:30,label:x_name,textsize:13},
                y:{size:45,label:y_name,textsize:13}
            }
        }
        super(dataStore,div,config,{x:{},y:{}});
        if (config.viv){
            this.addMenuIcon("fas fa-palette","Alter Channels")
            .addEventListener("click",(e)=>new ColorChannelDialog(this.viv,this.__doc__));
        }
       
    }

   
    setSize(x,y){
        super.setSize(x,y);
        if (this.viv){
            const b = this._getContentDimensions();
            this.viv.setSize(b.width,b.height,this.app);

        }
        

    }

    // PJT *Channel methods moved to VivViewer.

    remove(){
        this.viv.deck.finalize();
        super.remove();  
    }

    centerGraph(){
        super.centerGraph();
        if (this.viv){
            this.viv.setPanZoom(this.app.offset,this.app.x_scale,this.app.y_scale)
        }
    }

    getConfig(){
        const conf = super.getConfig();
        const k = this.viv.layers[0].props;
        conf.viv.image_properties={
            selections:k.selections.slice(0),
            colors:k.colors.slice(0),
            channelsVisible:k.channelsVisible.slice(0),
            contrastLimits:k.contrastLimits.slice(0),
            domains: k.domains?.slice(0)
        }
        return conf;
    }

    afterAppCreation(){
        const c = this.config;
        //make sure svg is on top of scatter plot and events pass through
        //this.contentDiv.prepend(this.app.div_container);
        // this.svg.style("position","absolute")
        //       .style("pointer-events","none");
        if (!c.viv) {
            // this dummy config isn't really adequate
            c.viv = {
                url: 'https://viv-demo.storage.googleapis.com/LuCa-7color_3x3component_data.ome.tif',
                image_properties: {
                    contrastLimits: [[0, 255]],
                    colors: [[255, 255, 255]],
                    channelsVisible: [true],
                    selections: [{c: 0, t: 0, z: 0}],
                }
            }
        }
        if (c.viv){ 
            const box =this._getContentDimensions();
        
            this.vivCanvas= createEl("canvas",{
                height:box.height,
                width:box.width,
                styles:{
                    position:"absolute",
                }
            });
            this.graphDiv.prepend(this.vivCanvas);
            

            import ('../webgl/VivViewer.js').then(({default:VivViewer})=>{
                this.viv = new VivViewer(this.vivCanvas,c.viv,this.app);
            
                this.app.addHandler("pan_or_zoom",(offset,x_scale,y_scale)=>{
                    this.viv.setPanZoom(offset,x_scale,y_scale)
                },c.id+"_viv")
            });
        }

        return super.afterAppCreation();
    }
}



BaseChart.types["viv_scatter_plot"]={
    name:"Viv Scatter Plot",
    allow_user_add:false,
    class:VivScatterPlot,
    params:[{
        type:"number",
        name:"X axis"
    },
    {
        type:"number",
        name:"Y axis"
    },
    {
        type: "text",
        name: "Category Column"
    }
    ]
}

export default VivScatterPlot;