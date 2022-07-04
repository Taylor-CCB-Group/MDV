
import { createEl,getElDim } from "./Elements";

class AutoComplete{
    constructor(element,getChoices,config={}){
        this.__doc__= config.doc || document;
        this.input=element;
        this.input.addEventListener("input", e => {
            this.text= this.input.value;
            getChoices(this.text).then(data=>{
                this.list=data;
                this.createDropDown();
            });      
        });
    }

    addItemEventListener(div, item){
        div.addEventListener("click", e=> {
           this.listHolder.remove();
           this.selectedItem=item;
           this.input.value=item.name      
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
                width:dim.width +"px"
            }
        },this.__doc__.body);
        
     
        /*for each item in the array...*/
        for (let i = 0; i < this.list.length; i++) {
            /*check if the item starts with the same letters as the text field value:*/
            const label = this.list[i].name;
              /*create a DIV element for each matching element:*/
              let b = createEl("div",{text:label},this.listHolder);
              this.addItemEventListener(b,this.list[i])
            
          }

    }
}

export default AutoComplete;