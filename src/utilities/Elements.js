
import Split from "split.js"

function createEl(type,attrs,parent){
   
    const el = document.createElement(type);
    
    if (attrs){
        addElProps(el,attrs)
    } 
    if (parent){
        parent.append(el);
    }
    return el;
}

function createSVGEl(type,attrs,parent){
  
    const el = document.createElementNS("http://www.w3.org/2000/svg", type);

    if (attrs){
        for (var idx in attrs) {
            if ((idx === 'styles' || idx === 'style') && typeof attrs[idx] === 'object') {
                for (var prop in attrs[idx]){el.style[prop] = attrs[idx][prop];}
            } else if (idx === 'text') {
                el.textContent = attrs[idx];
            } else if (idx==="classes"){
                for (var cl of attrs[idx]){el.classList.add(cl)}
                 
            } else {
                el.setAttributeNS(null,idx, attrs[idx]);
            }
        }
    } 
    if (parent){
        parent.append(el);
    }
    return el;
}

function splitPane(el,config={}){
    const dir = config.direction || "horizontal";
    const number = config.number || 2;
    const panes = [];
    const classes= dir==="horizontal"?["split-horizontal"]:["split-vertical"];
    for (let i =0;i<number;i++){
        panes.push(createEl("div",{classes:classes},el));
    }
    Split(panes,{
        direction:dir,
        gutterSize:5
    })
    return panes;
 
}

function createMenuIcon(icon,config,parent){
    const attrs={
    };
    const t  =config.tooltip;
    if(t){
        Object.assign(attrs,{
            "aria-label":t.text,
            role:"tooltip",
            "data-microtip-size":t.size || "small",
            "data-microtip-position":t.position || "bottom-left"
        });
    }

    const sp= createEl("span",attrs);

    createEl("i",{  
        classes:["ciview-menu-icon"].concat(icon.split(" ")),
        styles:{
            fontSize: config.size || "18px",
          
        }
    },sp);
    if (config.func){
        sp.addEventListener("click",(e)=>config.func(e));
    }
    if (parent){
        parent.append(sp);
    }
    return sp;

}

function addElProps(el,attrs){
    for (var idx in attrs) {
        if ((idx === 'styles' || idx === 'style') && typeof attrs[idx] === 'object') {
            for (var prop in attrs[idx]){el.style[prop] = attrs[idx][prop];}
        } else if (idx === 'text') {
            el.textContent = attrs[idx];
        } else if (idx==="classes"){
            for (var cl of attrs[idx]){el.classList.add(cl)}
             
        } else {
            el.setAttribute(idx, attrs[idx]);
        }
    }
}



function makeResizable(el,config={}){
    //already resizable
    if (el.__resizeinfo__){
        return;
    }
    const ri = {
        resize: el.style.resize,
        overflow:el.style.overflow
    }
    el.style.resize="both";
    el.style.overflow="hidden";
    //el.style.zIndex="0";
    if (config.onresizeend){
        ri.onresize=addResizeListener(el,(x,y)=>{
            config.onresizeend(x,y);
        })
    }
    el.__resizeinfo__=ri;
}

function removeResizable(el){
    if (!el.__resizeinfo__){
        return;
    }
    const ri = el.__resizeinfo__
    el.style.resize=ri.resize;
    el.style.overflow=ri.overflow;
    if (ri.onresize){
        el.removeEventListener("mouseup",ri.onresize)
    }
    delete el.__resizeinfo__;


}

function removeDraggable(el){
     const d = el.__draginfo__;
     if (!d){
         return;
     }
     d.handle.style.cursor=d.cursor;
     el.style.position=d.position;
     d.handle.onmousedown=null;
     delete el.__draginfo__;
}



function makeDraggable(el,config={}){
    if (!config.doc){
        config.doc=document;
    }
    let pos1 = 0, pos2 = 0, pos3 = 0, pos4 = 0;
    let  handle= config.handle?el.querySelector(config.handle):el;
    let cont = null;
    let is_moving=false;
    if (config.contain){
        cont={
            dir:config.contain
        };
    }



    handle.onmousedown = dragMouseDown;
    el.__draginfo__={
        handle:handle,
        cursor:handle.style.cursor,
        position:el.style.position,
    }
    handle.style.cursor="move";
    el.style.position="absolute";
    el.style.margin="0px 0px 0px 0px";
    el.__doc__=config.doc;
    function dragMouseDown(e) {
      e = e || window.event;
      e.preventDefault();
      // get the mouse cursor position at startup:
      pos3 = e.clientX;
      pos4 = e.clientY;
      if (config.ondragstart){
          config.ondragstart();
      }
      if (cont){
          cont.p_bb= el.parentElement.getBoundingClientRect();
          cont.c_bb =el.getBoundingClientRect();
      }
      el.__doc__.onmouseup = closeDragElement;
      el.__doc__.onmousemove = elementDrag;
    }

    
  
    function elementDrag(e) {
      e = e || window.event;
      e.preventDefault();
      e.stopPropagation();
      // calculate the new cursor position:
      pos1 = pos3 - e.clientX;
      pos2 = pos4 - e.clientY;
      pos3 = e.clientX;
      pos4 = e.clientY;
      //console.log(pos3+":"+pos4);
      console.log(pos1+":"+pos2);
      // set the element's new position:
      let nt  = el.offsetTop -pos2;
      let nl = el.offsetLeft - pos1;
    if (cont){
        if (nt<0 || (nt+cont.c_bb.height>cont.p_bb.height && cont.dir !=="topleft")){
            return;
        }
        if (nl<0 || (nl+cont.c_bb.width>cont.p_bb.width && cont.dir !=="topleft")){
            return;
        }
    }

      el.style.top = (nt) + "px";
      el.style.left = (nl) + "px";
    }
  
    function closeDragElement() {
      // stop moving when mouse button is released:
      if (config.ondragend){
          config.ondragend()
      }
      el.__doc__.onmouseup = null;
      el.__doc__.onmousemove = null;
    }
}

function addResizeListener(element,callback){
    let box = element.getBoundingClientRect();
    let width = box.width;
    let height = box.height;
    const list = (e)=>{
        let box = element.getBoundingClientRect();
        if (box.width!==width || box.height!==height){
            callback(box.width,box.height);
        }
        width=box.width;
        height=box.height;
    }
    element.addEventListener("mouseup",list)
    return list;
}

export {createEl,createSVGEl,addResizeListener,makeDraggable,
    makeResizable,removeDraggable,removeResizable,createMenuIcon,
    splitPane};