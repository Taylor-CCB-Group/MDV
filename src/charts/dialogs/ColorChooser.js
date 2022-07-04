import {BaseDialog} from "../../utilities/Dialog.js";
import {createEl} from "../../utilities/Elements.js";

class ColorChooser extends BaseDialog{
    constructor(cm){
        const config={
            footer:true,
            width:500,
            maxHeight:500,
            title:"Field Colors",
            columns:2,
            buttons:[{
                text:"Apply",
                method:"applyColor"
            }]
        }
        super(config,cm);
        
    }

    init(cm){
        this.cm =cm;
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
        createEl("label",{text:"Data Set:",styles:lStyles},d1)
        this.dsSelect= createEl("select",{
            classes:["mdv-select"]
        },d1);
        const dsNames= [];
        for (let d in cm.dsIndex){
            dsNames.push(d);
            createEl("option",{
                text:d,
                value:d
            },this.dsSelect)
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

        this.setDataSource(dsNames[0])
        createEl("span",{
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
        const cols = ds.getColumnList("text");
        for (let c of cols){
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
      const colors  = this.colorChoosers.map(x=>x.value);
      this.dataSource.setColumnColors(this.column,colors);
      this.dataSource.dataChanged([this.column],false,false);
      for (let ds of this.cm.dataSources){
        if (ds.column_link_to && ds.column_link_to.dataSource === this.dsName){
            this.cm._sync_colors(this.cm.dsIndex[ds.column_link_to.dataSource],ds);
        }
    }
    }

    setColumn(col){
        const c1= this.columns[0];
        this.column=col;
        c1.innerHTML="";
        const vals = this.dataSource.getColumnValues(col);
        const colors = this.dataSource.getColumnColors(col);
        this.colorChoosers=[];
        for (let n=0;n<colors.length;n++){
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
                text:vals[n]
            },d);
            const cc= createEl("input",{
                classes:["mdv-flex-fixed"],
                styles:{
                    width:"25px",
                    height:"20px",
                    padding:"0px",
                },
                type:"color",
                value:colors[n]
            },d);
            this.colorChoosers.push(cc);
        }
    }
}

export default ColorChooser;