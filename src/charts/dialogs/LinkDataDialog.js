import {BaseDialog} from "../../utilities/Dialog.js";
import {createEl} from "../../utilities/Elements.js";
import {getRandomString} from "../../utilities/Utilities.js";

class LinkDataDialog extends BaseDialog{
    constructor(cm,ds){
        const config={
            title:`Add Data To ${ds.link_to.dataSource}`,
            footer:true,
            width:380
        }
        super(config,{cm:cm,ds:ds});
    }
    init(content){
      
        this.cm=content.cm;
        this.ds=content.ds;
      
        this.link= this.ds.link_to;
        this.linkDataStore =  this.cm.dsIndex[this.link.dataSource].dataStore;
        this.dataTypes= this.linkDataStore.columnGroups[this.link.columnGroup].subgroups;
        this.ids=new Set();     
        createEl("div",
            {
                classes:["ciview-padded-div"],
                text:`Select row(s) in a table or click on point(s) in a plot to choose.
                Or paste in a list of ${this.link.choose_field.label}(s)`
        },this.dialog);

        const pf = createEl("div",{styles:{display:"flex",height:"150px",padding:"4px",gap:"8px"}},this.dialog);


        const gld= createEl("div",{styles:{flex:1,display:"flex",flexDirection:"column"}},pf);
        this.linputList= createEl("textarea",{},gld);
        //const b = createEl("button",{classes:["ciview-button-sm"],text:"ADD>"},gld);
       /* b.addEventListener("click",()=>{
            this.addToList();
        });*/
        const lld= createEl("div",{styles:{flex:1,display:"flex",flexDirection:"column"}},pf);
        createEl("div",{text:"Chosen"},lld);
        this.list = createEl("div",{styles:{flexGrow:1,overflowX:"hidden"},overFlowY:"scroll"},lld);
        this.listenerID=getRandomString();
        this.ds.dataStore.addListener(this.listenerID,(type,data)=>{
            if (type==="data_highlighted"){
                this.clearList();
                this.addIdsToList(data.indexes);
            }
            
        });

        const d= createEl("div",{classes:["ciview-padded-div"]},this.dialog);
        createEl("span",{classes:["ciview-title-div"],text:"Data Type:"},d);
        this.dataTypeSelect = createEl("select",{},d);

        //
        for (let id in this.dataTypes){
            const dt = this.dataTypes[id];
            createEl("option",{value:id ,text:dt.label},this.dataTypeSelect)
        }
       
        this._addColorChartSection();
        this._addHeatMapSection();

    }

    close(){
        this.ds.dataStore.removeListener(this.listenerID);
        super.close();
    }
  
      
        
    _addColorChartSection(){
        const div = createEl("div",{},this.dialog);
        createEl("div",{classes:["ciview-title-div"],text:"Color Chart"},div);
        
        createEl("span",{text:"Chart:",styles:{marginRight:"4px"}},div);
        const select = createEl("select",{styles:{maxWidth:"150px"}},div);
        const b= createEl("button",{classes:["ciview-button-sm"],text:"GO"},div);
        b.addEventListener("click",()=>{
            const param = this.getParams();
            if (param.length===0){
                return;
            }
            const ch = this.cm.getChart(select.value);
            ch.colorByColumn(param[0]);
        })

        for (let cid in this.cm.charts){
            const cinfo= this.cm.charts[cid];
            if (cinfo.dataSource.name===this.link.dataSource){
                const ch = cinfo.chart;
                if (ch.getColorOptions().colorby){
                    createEl("option",{value:cid,text:ch.config.title},select);
                }
            } 
        }

    }

    _addHeatMapSection(){
        const div = createEl("div",{},this.dialog);
        createEl("div",{classes:["ciview-title-div"],text:"Add Heat Map"},div);
        
        createEl("span",{text:"Category:",styles:{marginRight:"4px"}},div);
        const select = createEl("select",{styles:{maxWidth:"150px"}},div);
        const b= createEl("button",{classes:["ciview-button-sm"],text:"GO"},div);
        b.addEventListener("click",()=>{
            let param = this.getParams();
            if (param.length===0){
                return;
            }
            param= [select.value].concat(param);
            const dt = this.dataTypes[this.dataTypeSelect.value];
            this.cm.addChart(this.link.dataSource,{
                param:param,
                type:"heat_map",
                title:dt.label
            });
            
        });

        const list = this.linkDataStore.getColumnList("text");
        for (let item of list){
                
            createEl("option",{value:item.field,text:item.name},select);           
        }
    }


    getParams(){
        //check to see if already have metadada
        const dt = this.dataTypeSelect.value
        const fieldIds=[];
        for (let id of this.ids){
            const cf = this.link.choose_field.id;
            const name  = this.ds.dataStore.getRowText(id,cf);
            fieldIds.push(`${dt}|${name}(${dt})|${id}`);
           
        }
        return fieldIds;     
    }    
    
    clearList(){
        this.ids.clear();
        this.list.innerHTML="";
    }


    addIdsToList(ids){
        const f = this.link.choose_field.id;
        for (let id of ids){
            if (this.ids.has(id)){
                continue;
            }
            this.ids.add(id);
            const name = this.ds.dataStore.getRowAsObject(id,[f]);
            const item = createEl("div",{},this.list);
            const sp=createEl("span",{},item);
            sp.innerHTML=name[f];
        }
    }
}

export default LinkDataDialog;