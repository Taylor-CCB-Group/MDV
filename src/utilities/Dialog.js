import {createEl,makeResizable,makeDraggable} from "./Elements.js";

class BaseDialog{
  /**
  * Adds a menu icon with tooltip to the title bar 
  * @param {object} config The css classs of the icon (space delimited).
  * @param {string} [config.title] The dialog title
  * @param {DOMElement} [config.doc=document] The document to attach the dialog to 
  * @param {integer} [config.width] - The width of the dialog (in pixels)
  * @param {integer} [config.height] - The height of the dialog (in pixels)
  * @param {integer} [config.maxHeight] - if no height is specified,
  * the the dialog will fit its contents, unless maxHeight is given
  * @param {integer[]} [config.position] - The x and y psoition (in the document's body)
  * If absent, the dialog will be placed in the middle of the screen
  * @param {integer} [config.columns] An integer which will divide the dialog into
  * equally sized vertical columns, which can be accessed with this.columns array
  * @param {boolean} [config.footer=false] If true - a footer will be added and can be 
  * accessed with this.footer
  * @param {object[]} [config.buttons] -A list of of objects which should contain text
  * (the button's label) and method (the name of  the method to call when the button is clicked)
  *  e.g [{text:"OK",method:"doSomething"}] - the method 'doSomething' needs to be in the 
  * subclass. Buttons are added to the footer, so this needs to also be specified in the config
  * @param {function} [onClose] A function called when the dialog is closed
  * @param {object} content The object passed to the init method
  */

constructor (config={},content) {
    config.doc=config.doc || document;
    this.config=config;
    this.buttons={};
    const width = config.width?config.width+"px":"";
    const height = config.height?config.height+"px":"";
   

    this.outer= createEl("div",{
        classes:["ciview-dlg-outer"],
        styles:{
          width:width,
          height:height
        }
    });
    this.header = createEl("div",{
        classes:["ciview-dlg-header"]   
    },this.outer);

    createEl("span",{
        classes:["ciview-dlg-header-text"],
        text:config.title
    },this.header)

    createEl("i",{
        classes:["fas","fa-times","ciview-dlg-close-icon"]
    },this.header)
    .addEventListener("click",()=>{
      this.close();
    });

    this.dialog=createEl("div",{
      classes:["ciview-dlg-content"]
    },this.outer);

    makeResizable(this.outer,{
      doc:config.doc,
      onresizeend:(x,y)=>{
      this.onResize(x,y)
    }});
    makeDraggable(this.outer,{
      doc:config.doc,
      handle:".ciview-dlg-header"    
    });
   
    if (config.columns){
      this._createColumns(config.columns);
    }

    if (config.footer){
      this._addFooter();
    }
    if (config.buttons){
      if (!this.footer){
        this._addFooter();
      }
      for (let but of config.buttons){
        this._addButton(but);
      }
    }

    this.init(content);

    //may need to adjust its position depending on size
    //to avoid it being off screen
    config.doc.body.append(this.outer);
    const bbox = config.doc.body.getBoundingClientRect();
   
    const dbox =this.outer.getBoundingClientRect();
    if (config.maxHeight && dbox.height>config.maxHeight){
      dbox.height = config.maxHeight;
      this.outer.style.height=config.maxHeight+"px";
    }


    let pos = this.config.position;
    
    //put it in the middle of the screen if no position specified
    if (!pos){
      pos=[
        bbox.width/2 - dbox.width/2,
        bbox.height/2 - dbox.height/2
      ]
    }

   
    if (pos[0]+dbox.width > bbox.width){
      let w= bbox.width-dbox.width;
      pos[0]=w<0?0:w;
      
    }
    if (pos[1]+dbox.height > bbox.height){
      let h= bbox.height-dbox.height;
      pos[1]=h<0?0:h;
      
    }
    this.outer.style.left =pos[0]+"px";
    this.outer.style.top =pos[1]+"px"
    

  }
  setParent(parent) {
    if (parent) {
      parent.append(this.outer);
    } else {
      this.config.doc.body.append(this.outer);
    }
  }

  _addFooter(){
    this.footer= createEl("div",{
      styles:{
        flex:"0 0 auto",
        padding:"5px",
        justifyContent:"center",
        display:"flex",
        gap: "1em"
      }
    },this.outer);
  }

  _addButton(but){
    const b = createEl("button",{
      text:but.text,
      classes:["ciview-button"]
    },this.footer)
    b.addEventListener("click", ()=>this[but.method]());
    if (but.id){
      this.buttons[but.id]=b;
    }
  }

  disableButton(id,disable){
    //this.buttons[id].classList.remove("ciview-button");
    if (disable){
      this.buttons[id].classList.add("ciview-disabled");
    }
    else{
      this.buttons[id].classList.remove("ciview-disabled");
    }
    
  }

  _createColumns(number){
    this.dialog.style.display="flex";
    this.dialog.style.flexDirection="row";
    this.columns=[];
    for (let n=0;n<number;n++){
        this.columns.push(createEl("div",{
          styles:{
            height:"100%",
            flexBasis:"0",
            flexGrow:"1",
            padding:"0px 4px",
            flexDirection:"column",
            alignItems:"stretch",
          }
        },this.dialog));
    }
  }

  _closed = false;
  /**
  * Closes the dialog and removes it from the DOM
  */
  close(){
    if (this._closed) return;
    this._closed = true;
    if (this.config.onclose){
      this.config.onclose();
    }
    this.outer.remove()

  }

  /**
  * This method should be overridden to add content to the dialog
  */
  init(content){}

  /**
  * This method should be overridden if the layout needs
  * updating on resize
  * @param {integer} x the new width of the dialog
  * @param {integer} y the new height of the dialog
  */
  onResize(x, y) {
    if (this.outer) {
        this.outer.style.width = x + "px";
        this.outer.style.height = y + "px";
    }
}


}




function getTextInput(title,event,doc=document){
  return new Promise((accept,reject)=>{
    const div= createEl("div",{
      styles:{
        position:"absolute",
        zIndex:200,
        left:event.clientX+"px",
        top:event.clientY+"px",
        background:"white",
        padding:"5px"
      }
  
    },doc.body)
    createEl("label",{text:title,style:{display:"block"}},div);
    const i = createEl("input",{},div);
    i.addEventListener("keypress",(e)=>{
      if (e.keyCode===13){
       submit();
      }
    });
    i.focus();
    const el= ()=>submit();
    i.addEventListener("blur",el);

    function submit(){
      i.removeEventListener("blur",el);
      div.remove();
      accept(i.value);
      
    }

  })


}

//move this, type it differently... figure out what it is that makes it tick, or not.
BaseDialog.experiment = {};

export {BaseDialog,getTextInput};