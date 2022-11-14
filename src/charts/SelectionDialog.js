import BaseChart from "./BaseChart.js"
import { createEl } from "../utilities/Elements.js";
import noUiSlider from "nouislider";

class SelectionDialog extends BaseChart{
    constructor(dataStore,div,config){
        super(dataStore,div,config);
        const c  =this.config;
        const con= c.param;
        this.dims={};
        this.textFilters={};
        this.numberFilters={};
        c.filters=c.filters || {};
        this.hasFiltred=false;
        for (let c of con){
            const col =dataStore.columnIndex[c];
            const div =createEl("div",{
                styles:{padding:"5px"}
            },this.contentDiv)
            createEl("div",{text:col.name,styles:{fontWeight:"bold"}},div)
            if (col.datatype==="text"){
                this.addTextFilter(col,div)
            }
            else if (col.datatype==="integer" || col.datatype==="double"){
                this.addNumberFilter(col,div);
            }
        }

        if (this.hasFiltered){
            setTimeout(()=>this.dataStore.triggerFilter(),1500);
        }
    }

    removeFilter(){
        const f= this.config.filters   
        for (let c in this.textFilters){
            const t = this.textFilters[c];
            t[0].value="__none__";
            t[1].checked=false;
            f[c].category="__none__";
            f[c].excluse=false;
        }
        for (let c in this.numberFilters){
            const mm= this.dataStore.getMinMaxForColumn(c);
            const t = this.numberFilters[c];
            t.noUiSlider.set([mm[0],mm[1]])
            f[c]=null;
        }
        this.resetButton.style.display="none";
        for (let c in this.dims){
            this.dims[c].removeFilter(false);
        }
        this.dataStore.triggerFilter();
        
    }

    onDataFiltered(dim){
     
       
    }

    addNumberFilter(col,div){
        const sl = createEl("div",{
            styles:{
                margin:"0px 10px"
            }
        },div);
        const c = this.config;     
        const dim =this.dataStore.getDimension("range_dimension");
        dim.noclear=true;
        this.dims[col.field]= dim;
        const mm = this.dataStore.getMinMaxForColumn(col.field)

        
        const fil = c.filters[col.field];
        let cv = [mm[0],mm[1]]
        if (fil){
            cv=[fil[0],fil[1]]
        }
        noUiSlider.create(sl, {
            start: [cv[0],cv[1]],
            range: {
                'min': mm[0],
                'max': mm[1]
            },
            //step:s.step || null,
            tooltips:true,
            documentElement:this.__doc__
        });
        sl.noUiSlider.on("end",values=>{
            const min = parseFloat(values[0]);
            const max = parseFloat(values[1]);
            if (min<=mm[0] && max >=mm[1]){
                c.filters[col.field]=null;
            }
            else{
                c.filters[col.field]=[min,max];
            }
            
            this.filterCategories(col.field);
        });
        if (fil){
            this.filterCategories(col.field,false);
            this.hasFiltered=true;
        }
        this.numberFilters[col.field]=sl;
    }

    filterCategories(col,notify=true){
        const f= this.config.filters[col];
        if (this.dataStore.columnIndex[col].datatype=== "text"){
            if (f.category === "__none__"){
                this.dims[col].removeFilter();
            }
            else{
                let vs= [f.category]
                if (f.exclude){
                    vs = this.dataStore.getColumnValues(col).filter(x=>x!==f.category)
                }
                this.resetButton.style.display = "inline";
                this.dims[col].filter("filterCategories",[col],vs,notify)
            }    
        }
        else{
            if (!f){
                this.dims[col].removeFilter();
            }
            else{
                this.resetButton.style.display = "inline";
                this.dims[col].filter("filterRange",[col],{min:f[0],max:f[1]},notify)
            }

        }
        
    }

    addTextFilter(col,div){
        div.style.whiteSpace="nowrap";
        const c = this.config;     
        const dim =this.dataStore.getDimension("category_dimension");
        dim.noclear=true;
        this.dims[col.field]= dim;

        if (!c.filters[col.field]){
            c.filters[col.field]={
                exclude:false,
                category:"__none__"
            }
        }
        else{
            this.filterCategories(col.field,false);
            this.hasFiltered=true;
        }
        const fil = c.filters[col.field]
        const dd = createEl("select",{
            classes:["mdv-select"],
            styles:{
                maxWidth:"200px"
            }
        });
        createEl("option",{
            text:"None",
            value:"__none__"
        },dd)
        for (let name of col.values){
            createEl("option",{
                text:name,
                value:name
            },dd)
        }
        div.append(dd);
        dd.value = fil.category
        dd.addEventListener("change",(e)=>{
            fil.category=dd.value;
            this.filterCategories(col.field);
          
        });
        const cb =createEl("input",{
            classes:["mdv-checkbox"],
            type:"checkbox"
        },div)
        cb.checked= fil.exclude;

        cb.addEventListener("click",e=>{
            fil.exclude = cb.checked;
            this.filterCategories(col.field);     
        });
        
        createEl("span",{
            text:"exclude",
            style:{
                fontSize:"11px",
                verticalAlign:"middle"
            }
        },div);
        this.textFilters[col.field]=[dd,cb]   
        
    }


    remove(notify=true){
        for (let c in this.dims){
            this.dims[c].destroy(false);
        }
        if (notify){
            this.dataStore.triggerFilter();
        }
        super.remove();
    }
}

BaseChart.types["selection_dialog"]={
    name:"Selection Dialog",
    "class":SelectionDialog,
    params:[
        {
        type:"_multi_column:all",
        name:"Columns To filter"
        }
    ]
}

export default SelectionDialog;