import {BaseDialog} from "../../utilities/Dialog.js";
import {createEl} from "../../utilities/Elements.js";

class ColorChooser extends BaseDialog{
    constructor(cm,ds){
        const config={
            footer:true,
            width:500,
            maxHeight:500,
            title:"Field Colors For "+ds.name,
            columns:2,
            buttons:[{
                text:"Apply",
                method:"applyColor"
            }]
        }
        super(config,{cm:cm,ds:ds});
        
    }

    init(content){
        const cm = this.cm =content.cm;
        this.ds = content.ds;
        const c2 = this.columns[1];
        this.columns[0].style.overflowY="auto";
        const dStyle={
            whiteSpace:"nowrap",
            padding:"2px"

        }
        const d1 = createEl("div",{styles:dStyle},c2);
        const lStyles={
            width:"70px",
            display:"inline-block"
        }
       
        const d2 = createEl("div",{styles:dStyle},c2);
        createEl("label",{text:"Field:",styles:lStyles},d2)
        this.columnSelect = createEl("select",{
            classes:["mdv-select"],
            styles:{
                maxWidth:"150px"
            }
        },d2);

        const d3 = createEl("div",{styles:dStyle},c2);
        createEl("div",{
            text:"Color Scheme",
            styles:{
                fontWeight:"bold"
            }
        },d3)
        createEl("span",{
            classes:["mdv-info-text"],
            text:"Paste list of colors (hex)"
            
        },d3);

        const d4 = createEl("div",{
            styles:{
                display: "flex",
                alignItems: "flex-start",
                marginTop:"2px"
            }

        },c2)
        this.colorList=createEl("textarea",{
            classes:["mdv-select"],
            styles:{
                width:"100px",
                padding:"2px",
                height:"150px"
            }
        },d4);

        this.setDataSource(this.ds.dataStore.name);
        createEl("button",{
            classes:["ciview-button-sm"],
            text:"Use Scheme"
        },d4).addEventListener("click",()=>{
            this.parseColors(this.colorList.value);
        })
    }
    parseColors(list){
        const hex =/^#[0-9A-F]{6}$/i;

        let colors  = list.replaceAll("\"","").split(/[,\s]+/);
        colors=colors.filter(x=>x.match(hex))
        const colorList= colors.join("\n");
        this.colorList.value=colorList;
        let index = 0
        if (colorList.length>0){
            for (let n=0;n<this.colorChoosers.length;n++){
                this.colorChoosers[n].value= colors[index];
                index++;
                if (index===colors.length){
                    index=0;
                }
            }
        }

    }

    setDataSource(name){
        this.columnSelect.innerHTML ="";
        const ds = this.cm.getDataSource(name);
        this.dataSource=ds;
        this.dsName=name;
        const cols = ds.getColumnList("string");
        for (const c of cols){
            createEl("option",{
                text:c.name,
                value:c.field
            },this.columnSelect)
        }

        this.columnSelect.addEventListener("change",()=>{
            this.setColumn(this.columnSelect.value);
        });
        this.setColumn(cols[0].field)

    }

    applyColor(){   
        const colors = new Array(this.colorChoosers.length)
        this.colorChoosers.forEach(x=>colors[x.__index]=x.value);
        this.dataSource.setColumnColors(this.column,colors);
        this.dataSource.dataChanged([this.column],false,false);
        for (const d of this.cm.dataSources){
            const ds  = d.dataStore;
            const scc = ds.syncColumnColors.find(x=>x.dataSource===this.dataSource.name);
            if (scc){
                const  c= scc.columns.find(x=>x.link_to===this.column);
                if (c){
                    this.cm._sync_colors([c],this.dataSource,ds);
                    ds.dataChanged([c.col],false,false);     
                }
            }
            const lcc = ds.linkColumns.find(x=>x.dataSource===this.dataSource.name);
            if (lcc){
                if (lcc.columns.indexOf(this.column)!==-1){
                    ds.setColumnColors(this.column,colors)
                    ds.dataChanged([this.column],false,false);            
                }      
            }
        }
    }

    setColumn(col){
        const c1= this.columns[0];
        this.column=col;
        c1.innerHTML="";
        const vals = this.dataSource.getColumnValues(col);
        const sorted = vals.slice(0).sort();
        const colors = this.dataSource.getColumnColors(col);
        const sortedColors =  colors.map((x,i)=>{
            const index = vals.indexOf(sorted[i]);
            return {index:index,color:colors[index]}
        })
        this.colorChoosers=[];
        for (let n=0;n<sorted.length;n++){
            const d = createEl("div",{styles:{
                display:"flex",
                flexWrap:"nowrap"
            }},c1);
            createEl("span",{
                styles:{
                    display:"inline-block",
                    fontSize:"0.8erm"
                },
                classes:["mdv-flex-dynamic"],
                text:sorted[n]
            },d);
            const cc= createEl("input",{
                classes:["mdv-flex-fixed"],
                styles:{
                    width:"25px",
                    height:"20px",
                    padding:"0px",
                },
                type:"color",
                value:sortedColors[n].color
            },d);
            cc.__index=sortedColors[n].index;
            this.colorChoosers.push(cc);
        }
    }
}

export default ColorChooser;