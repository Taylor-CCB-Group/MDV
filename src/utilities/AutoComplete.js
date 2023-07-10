
import { createEl,getElDim } from "./Elements";

class AutoComplete{
    constructor(element,getChoices,config={}){
        this.__doc__= config.doc || document;
        this.input=element;
        this.listeners=[];
        config.minLength= config.minLength || 2;
        this.input.addEventListener("input", e => {
            if (this.input.value.length<config.minLength){
                return;
            }
            this.text= this.input.value;
            getChoices(this.text).then(data=>{
                this.list=data;
                this.createDropDown();
            });      
        });

        //click outside dropdown removes the list

        this.dli = e=>{
            if (e.target.parentElement=== this.listHolder){
                return;
            }
            if (this.listHolder ){
                this.listHolder.remove();
                this.listHolder=null;
            }
        }
        this.__doc__.addEventListener("pointerdown",this.dli);
    }

    getSelectedItem(){
        return this.selectedItem;
    }

    addListener(func){
        this.listeners.push(func)
    }

    destroy(){
        if (this.listholder){
            this.listholder.remove();
        }
        this.__doc__.removeEventListener("click",this.dli);
    }

    setSelectedItem(item,notify=true){
        this.selectedItem=item;
        this.input.value=item.value;
        if (notify){
            this.listeners.map(x=>x("item_chose",item))
        }
    }

    addItemEventListener(div, item){
        div.addEventListener("click", e=> {
           this.listHolder.remove();
           this.listHolder=null;
           this.setSelectedItem(item);
        });
    }

    createDropDown(){
        if (this.listHolder){
            this.listHolder.remove();
        }
        if (!this.text){
            return;
        }
        const dim = getElDim(this.input);
        this.listHolder = createEl("div",{
            classes:["autocomplete-items"],
            styles:{
                top:(dim.top+dim.height)+"px",
                left:dim.left+"px",
                width:dim.width +"px",
                zIndex:1200
            }
        },this.__doc__.body);
        
   
        for (let i = 0; i < this.list.length; i++) {
            const label = this.list[i].value;
            let b = createEl("div",{text:label},this.listHolder);
            this.addItemEventListener(b,this.list[i])  
          }
    }
}

export default AutoComplete;