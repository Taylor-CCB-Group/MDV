import {createEl,makeResizable,makeDraggable} from "./Elements.js";

class BaseDialog{

constructor (config={},content) {
    config.doc=config.doc?config.doc:document;
    this.config=config;
    const width = config.width?config.width+"px":"";
    const height = config.height?config.height+"px":"";
   

    this.outer= createEl("div",{
        classes:["ciview-dlg-outer"],
        styles:{
          width:width,
          height:height
        }
    })
    let self = this;
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

  },this.outer)

    makeResizable(this.outer,{
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

    this.init(content);

    //may need to adjust its position depending on size
    //to avoid it being off screen
    config.doc.body.append(this.outer);
    const bbox = config.doc.body.getBoundingClientRect();
   
    const dbox =this.outer.getBoundingClientRect();
    if (config.maxHeight && dbox.height>config.maxHeight){
      dbox.height = config.MaxHeight;
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

  _addFooter(){
    this.footer= createEl("div",{
      styles:{
        flex:"0 0 auto",
        padding:"5px",
        justifyContent:"center",
        display:"flex"

      }
    },this.outer);
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
            padding:"0px 4px"
          }
        },this.dialog));
    }


  }

  close(){
    
    if (this.config.onclose){
      this.config.onclose();
    }
    this.outer.remove()

  }


  init(){}

  onResize(x,y){}


}

export {BaseDialog};