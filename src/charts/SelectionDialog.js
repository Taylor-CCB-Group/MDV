import BaseChart from "./BaseChart.js"
import { createEl,makeSortable} from "../utilities/Elements.js";
import noUiSlider from "nouislider";
import { getRandomString } from "../utilities/Utilities.js";


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
        this.contentDiv.style.overflowY="auto"
        for (let c of con){
            const col =dataStore.columnIndex[c];
            const div =createEl("div",{
                styles:{padding:"5px"}
            },this.contentDiv);
            //marker to identify which div associated with which column
            div.__fcol__=col.field;
            createEl("div",{text:col.name,classes:["i-title"],styles:{fontWeight:"bold"}},div)
            if (col.datatype==="text" ){
                this.addTextFilter(col,div)
            }
            else if (col.datatype==="multitext"){
                this.addMultiTextFilter(col,div);
            }
            else if (col.datatype==="integer" || col.datatype==="double"){
                this.addNumberFilter(col,div);
            }
        }
        makeSortable(this.contentDiv,
            {
                handle:"i-title",
                //re-arrange params in config 
                sortEnded:(li)=>{        
                    this.config.param=li.map(x=>x.__fcol__);
                }
            });


        if (this.hasFiltered){
            setTimeout(()=>this.dataStore.triggerFilter(),0);
        }
    }

    removeFilter(){
        const f= this.config.filters;
      
        for (let c in this.textFilters){
            const t = this.textFilters[c];
            const is_multi = f[c].operand;
            if (!is_multi){
                t[0].value="__none__";
                t[1].checked=false;
            }
            else{
                t[1].innerHTML="";
            }
            
            f[c].category=is_multi?[]:"__none__";
            f[c].exclude=false;
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
                margin:"10px 10px"
            }
        },div);
        const dd = createEl("div",{styles:{
            padding:"3px 10px",
            overflow:"hidden"
        }},div)
        createEl("span",{text:">",styles:{"float":"left","fontWeight":"bold"}},dd);
        const greaterThan = createEl("input",{
            styles:{
                width:"70px",
                float:"left"
            }
        },dd);
        createEl("span",{text:"<",styles:{"float":"right","fontWeight":"bold"}},dd);
        const lessThan = createEl("input",{
            styles:{
                width:"70px",
                float:"right"
            }
        },dd);
      

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
        sl.noUiSlider.on("set",values=>{
            const min = parseFloat(values[0]);
            const max = parseFloat(values[1]);
            if (min<=mm[0] && max >=mm[1]){
                c.filters[col.field]=null;
            }
            else{
                c.filters[col.field]=[min,max];
            }
            lessThan.value=max;
            greaterThan.value= min;

            
            this.filterCategories(col.field);
        });
        
      
        lessThan.value=cv[1];
        greaterThan.value=cv[0];

        greaterThan.addEventListener("blur",e=>{
                let n = parseFloat(greaterThan.value);
               
                const cvi = sl.noUiSlider.get();
                if (isNaN(n) || n > cvi[1]){
                    return;
                }
                n  = n < mm[0]?mm[0]:n;
                sl.noUiSlider.set([n, cvi[1]],true,true);
        });
        greaterThan.addEventListener("keypress",e=>{
            if (e.key=== "Enter"){
                greaterThan.blur();
            }
        });

        lessThan.addEventListener("blur",e=>{
            
                let n = parseFloat(lessThan.value);
               
                const cvi = sl.noUiSlider.get();
                if (isNaN(n) || n < cvi[0]){
                    return;
                }
                n= n>mm[1]?mm[1]:n;
                sl.noUiSlider.set([cvi[0], n],true,true);
              
        });
        lessThan.addEventListener("keypress",e=>{
            if (e.key=== "Enter"){
                lessThan.blur();
            }
        });


        if (fil){
            this.filterCategories(col.field,false);
            this.hasFiltered=true;
        }
        this.numberFilters[col.field]=sl;
    }

    

    filterCategories(col,notify=true){
        const f= this.config.filters[col];
        const type = this.dataStore.columnIndex[col].datatype;
        if (type=== "text" ){
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
        else if (type==="multitext"){
            if (f.category.length===0){
                this.dims[col].removeFilter();
            }
            else{
                this.resetButton.style.display = "inline";
                const  args= f.category.slice(0);
                args.operand= f.operand;
                this.dims[col].filter("filterCategories",[col],args,notify)
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

    createTag(div,col,value,add=true){
        if (add){
            this.config.filters[col.field].category.push(value);
        }
        const d=createEl("div",{},div);
        createEl("span",{text:value},d);
        const i = createEl("i",{classes:["fas","fa-times"]},d);
        i.addEventListener("click",e=>{
            this.config.filters[col.field].category=this.config.filters[col.field].category.filter(x=>x!==value);
            console.log(this.config.filters[col.field].category);
            d.remove();
            this.filterCategories(col.field);  
        })
    }

    addMultiTextFilter(col,hdiv){
        let div = createEl("div",{style:{whiteSpace:"nowrap"}},hdiv)
        const c = this.config;     
        const dim =this.dataStore.getDimension("category_dimension");
        dim.noclear=true;
        this.dims[col.field]= dim;

        if (!c.filters[col.field]){
            c.filters[col.field]={
                operand:"or",
                category:[]
            }
        }
        else{
            this.filterCategories(col.field,false);
            this.hasFiltered=true;
        }

        const fil = c.filters[col.field];
        const dd = createEl("select",{
            classes:["mdv-select"],
            styles:{
                maxWidth:"200px"
            }
        });
       
        for (let name of col.values.slice(0).sort()){
            createEl("option",{
                text:name,
                value:name
            },dd)
        }
        div.append(dd);
        const b = createEl("button",{
            classes:["ciview-button-sm"],
            text:"Add"
        },div);
        const tdiv = createEl("div",{},div);
        b.addEventListener("click",e=>{
            if (fil.category.indexOf(dd.value)===-1){
                this.createTag(tdiv,col,dd.value);
                this.filterCategories(col.field);  
            }
            
        });
        div= createEl("div",{},hdiv)
        const rname = getRandomString();
        createEl("span",{text:"and"},div)
        const rb = createEl("input",{   
            type:"radio",
            name:rname,
            value:"and"
        },div);
        rb.checked = fil.operand==="and";
        rb.addEventListener("click",(e)=>{
            fil.operand="and";
            this.filterCategories(col.field);  

        });
        createEl("span",{text:"or"},div)
        const ra = createEl("input",{   
            type:"radio",
            name:rname,
            value:"or"
        },div);
        ra.checked = fil.operand==="or";
        ra.addEventListener("click",(e)=>{
            fil.operand="or";
            this.filterCategories(col.field);  

        })
        for (let c of fil.category){
            this.createTag(tdiv,col,c,false);
        }
      
        this.textFilters[col.field]=[dd,tdiv];
        
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
        const fil = c.filters[col.field];
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