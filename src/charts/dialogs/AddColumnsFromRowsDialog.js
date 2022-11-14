import {BaseDialog} from "../../utilities/Dialog.js";
import {createEl} from "../../utilities/Elements.js";
import AutoComplete from "../../utilities/AutoComplete.js";
import {getRandomString} from "../../utilities/Utilities.js";

const max_genes=1000;

class AddColumnsFromRowsDialog extends BaseDialog{
    constructor(ds,ds_to,link,cm){
        const config={
            title:`Add ${link.name}`,
            footer:true,
            columns:2,
            width:380
        }
        super(config,{ds,ds_to,link,cm});
    }
    init(content){
        this.ds=content.ds;
        this.ds_to=content.ds_to;
        this.link= content.link;
        this.cm=content.cm;

        //add single gene
        const h = createEl("div",{classes:["mdv-section"]},this.columns[0]);
        const ni = createEl("input",{},h);
        this.ac= new AutoComplete(ni,async (text)=>{
            return this.ds_to.dataStore.getSuggestedValues(text,this.link.name_column,8)
        });

        this.ac.addListener((type,data)=>this.addColumn(data));

        
        //add the radio buttons for type
        this.rn = getRandomString();
        for (let sg in this.link.subgroups){
            createEl("label",{
                text:this.link.subgroups[sg].label
            },this.footer);
            createEl("input",{
                type:"radio",
                name:this.rn,
                value:sg,
                checked:true
            },this.footer);
        }


        
        this._addOptionsSection();

        const b = createEl("div",{classes:["ciview-button-sm"],styles:[],text:"Apply"},this.columns[0]);
        b.addEventListener("click",()=>{
            this.addColumn(this.ac.getSelectedItem());
        });

        /*this.ta  = createEl("textarea",{},this.dialog);
        const bu = createEl("span",{classes:["ciview-button-sm"],text:"Add All"},this.dialog);
        bu.addEventListener("click",e=>{
            this.parseNames(this.ta.value)
        });*/

        this.cm.addListener(this.rn,(type,c,data)=>{
            switch(type){
                case "chart_added":
                    this._addColorChartOption(data);
                    break;
                case "chart_removed":
                    this._addColorChartOption(data,true);
                    break;
                case "filtered":
                    console.log("ss");
                    break;
                case "view_laoded":
                    this.close();
            }
        });

        this.itemList=[];
        const tds=  this.ds_to.dataStore;
        tds.addListener(this.rn,(type,data)=>{
            if (type==="data_highlighted"){
                const index = data.indexes[0];
                const value= tds.getRowText(index,this.link.name_column);
                this.ac.setSelectedItem({value,index})
            }
            else if (type==="filtered"){
                if (tds.filterSize<max_genes && tds.isFiltered()){
                    //this.addFilteredItems();
                }
            }
        });
    }

    addFilteredItems(){
        const tds=  this.ds_to.dataStore;
        for (let n=0;n<tds.size;n++){
            if (tds.filterArray[n]===0){
                this.itemList.push(n);
                createEl("div",{text:tds.getRowText(n,this.link.name_column)},this.columns[1]);
            }
        }
    }

    close(){
        this.cm.removeListener(this.rn);
        this.ds_to.dataStore.removeListener(this.rn);
        this.ac.destroy();
        super.close();
    }

    _addColorChartOption(chart,remove=false){
        const c = chart.config;
        if (chart.dataStore !== this.ds.dataStore || !(chart.getColorOptions().colorby)){
            return;
        }
        if (remove){
            const opt = this.colorChartSelect.querySelector(`[value='${c.id}']`);
            if (opt){
                opt.remove();
            }
            if (this.colorChartSelect.children.length===0){
                this.colorDiv.classList.add("mdv-disabled");
            }
            return;
        }
        createEl("option",{value:c.id,text:c.title},this.colorChartSelect);
        this.colorDiv.classList.remove("mdv-disabled");  
    }

    _addOptionsSection(){
        this.colorDiv = createEl("div",{classes:["mdv-section","mdv-disabled"]},this.columns[0]);  
        this.colorCheck= createEl("input",{type:"checkbox"},this.colorDiv);
        createEl("span",{text:"Color Chart:",styles:{marginRight:"4px"}},this.colorDiv);
        createEl("br",{},this.colorDiv);
        this.colorChartSelect = createEl("select",{styles:{maxWidth:"150px"}},this.colorDiv);
      
        for (let cid in this.cm.charts){
            this._addColorChartOption(this.cm.charts[cid].chart);
           
        }
        const d = createEl("div",{classes:["mdv-section"]},this.columns[0]);  
        this.histoCheck= createEl("input",{type:"checkbox"},d);
        createEl("span",{text:"Create Histogram:",styles:{marginRight:"4px"}},d);
    }

    

    parseNames(text){
        let names =  text.replaceAll("\"","").split(/[,\s]+/);
        names=names.filter(x=>x!=="");
        this.ta.value="";
        let g="";
        for (let n of names){
            const resp = this.ds_to.dataStore.getNearestMatch(n,this.link.name_column,1,1);
            if (resp.length>0){
                if (resp[0].mismatches===0){
                    g+=resp[0].value+"\n"
                }
                else{
                    g+=n+" - "+resp[0].value+"\n";
                }        
            }
            else{
                g+=n+"??\n"
            }
            
        }
        this.ta.value=g;
    }

    addColumn(c){
        const checked = document.querySelector(`input[name='${this.rn}']:checked`);
        const sg = checked.value;
        const f = `${sg}|${c.value}(${sg})|${c.index}`
    
        this.ds.dataStore.addColumnFromField(f);
        if (this.colorCheck.checked){
            const v= this.colorChartSelect.value;
            if (v){
                this.cm.getChart(v).colorByColumn(f);
                
            }
        }

        if (this.histoCheck.checked){
            this.cm.addChart(this.ds.name,{
                type:"bar_chart",
                param:f
            })
        }
    }

}

export default AddColumnsFromRowsDialog;