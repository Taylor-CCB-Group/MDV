
import Split from "split.js"

/**@template {keyof HTMLElementTagNameMap} T
 * @param {T} type
 * @param {{styles?:, classes?:string[], text?:string, [key: string]:string|string[]}=} attrs
 * @param {HTMLElement=} parent
 * @returns {HTMLElementTagNameMap[T]}
 */
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
        gutterSize:5,
        cursor:"pointer",
        sizes:config.sizes,
        onDragEnd:()=>window.dispatchEvent( new Event('resize') )
    })
    return panes;
 
}

function createMenuIcon(icon,config,parent){
    const attrs={
        role: "button", //a11y - we could use an actual button rather than span here, would need styling.
    };
    const t  =config.tooltip;
    if(t){
        Object.assign(attrs,{
            "aria-label":t.text,
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

/**
 * Enables the children of an HTML element to be sorted by the user
 * via drag and drop
 * @param {HTMLElement} list - an HTML element whose immediate children will
 * be capable of being sorted.
 * @param {object} config 
 * @param {string} [config.handle] - if each child element is to be only 
 * dragged by a handle, then this should be a class name that the handle 
 * contains.
 * @param {function} [config.sortEnded] - a function that is called when
 * a sorting action has finished. The function is supplied with the
 * current sorted list of child elements.
 */
function makeSortable(list,config={}){
    let dragged = null;
    let allowDrag=false;

    function handleDragOver(e){
        //place the dragged item before the current item
        this.after(dragged);
        e.preventDefault();
        return false;
    }
    
    function handleDragEnd(e){
        allowDrag=false;
        dragged.classList.remove("mdv-element-dragged");
        if (config.sortEnded){
            config.sortEnded(Array.from(list.children).slice(1));
        }
    }

    function handleMouseDown(e){
        //check whether item can be sorted
        if (!config.handle){
            allowDrag=true
        }
        else{
            if(e.target.classList.contains(config.handle)){
                allowDrag=true
            }
            else{
                allowDrag=false;
            }
       }  
    }

    function handleDragStart(e){
        if (!allowDrag){
            e.preventDefault();
            return false;
        }
        dragged =this;
        this.classList.add("mdv-element-dragged");
        //blank drag image as the swapping of items
        //will give the appearance of dragging
        e.dataTransfer.setDragImage(createEl("div",{}),0,0);
    }

    //add dummy element, allows items to be dragged to the beginning
    //of the list i.e. after the dummy element
    let firstDummyEl= createEl("div",  {styles:{height:"5px"}});
    list.prepend(firstDummyEl);
    firstDummyEl.addEventListener("dragover",handleDragOver);
    //add drag events to all child elements
    for (let n=0;n<list.children.length;n++){
        const item = list.children[n];
        item.draggable=true;
        item.addEventListener("dragover",handleDragOver);
        item.addEventListener("dragstart",handleDragStart);
        item.addEventListener("dragend",handleDragEnd);
        item.addEventListener("pointerdown",handleMouseDown);     
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
    // el.style.resize="both"; //standard resizer is sometimes visible when it shouldn't be.
    el.style.overflow="hidden";
    //el.style.zIndex="0";
    if (config.onresizeend){
        ri.onresize=addResizeListener(el,(x,y)=>{
            config.onresizeend(x,y);
        }, config.onResizeStart);
    }
    //workaround for Safari bug #50
    //https://codepen.io/jkasun/pen/QrLjXP
    //(actually, mostly copilot filling in very similar code...)
    const bottomRight = createEl("div", {
        classes: ["resizer-both"],
    }, el);
    bottomRight.addEventListener("pointerdown", initDrag, false);
    function initDrag(e) {
        ri.startX = e.clientX;
        ri.startY = e.clientY;
        ri.startWidth = parseInt(document.defaultView.getComputedStyle(el).width, 10);
        ri.startHeight = parseInt(document.defaultView.getComputedStyle(el).height, 10);
        document.documentElement.addEventListener("pointermove", doDrag, false);
        document.documentElement.addEventListener("pointerup", stopDrag, false);
    }
    function doDrag(e) {
        el.style.width = (ri.startWidth + e.clientX - ri.startX) + "px";
        el.style.height = (ri.startHeight + e.clientY - ri.startY) + "px";
    }
    function stopDrag(e) {
        document.documentElement.removeEventListener("pointermove", doDrag, false);
        document.documentElement.removeEventListener("pointerup", stopDrag, false);
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
        el.removeEventListener("pointerup",ri.onresize)
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


class MDVProgress{
    constructor(el,config){
        this.el=el;
        this.max=config.max || 100;
        el.classList.add("mdv-progress-outer");
        this.bar = createEl("div",{
            classes:["mdv-progress-inner"],
        },el);
        this.text=createEl("div",{
            classes:["mdv-progress-text"]
        },el);
        this.setText(config.text || "");
        this.setValue(config.value || 0 );
       
    }

    setText(text){
        this.text.textContent=text;  
    }
    setValue(val){   
        this.bar.style.width=`${Math.round((val/this.max)*100)}%`
    }
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
          config.ondragstart(el);
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
      //console.log(pos1+":"+pos2);
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
    if (!config.y_axis){
      el.style.top = (nt) + "px";
    }
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

function addResizeListener(element, endCallback, startCallback){
    let box = element.getBoundingClientRect();
    let width = box.width;
    let height = box.height;
    const list = (e)=>{
        let box = element.getBoundingClientRect();
        if (box.width!==width || box.height!==height){
            endCallback(box.width,box.height);
        }
        width=box.width;
        height=box.height;
    }
    if (startCallback) element.addEventListener("pointerdown", startCallback);
    element.addEventListener("pointerup",list)
    return list;
}

function getElDim(el){
    var rect = el.getBoundingClientRect(),
    scrollLeft = window.pageXOffset || document.documentElement.scrollLeft,
    scrollTop = window.pageYOffset || document.documentElement.scrollTop;
    return { top: rect.top + scrollTop, left: rect.left + scrollLeft ,height:rect.height,width:rect.width}

}

export {createEl,createSVGEl,addResizeListener,makeDraggable,
    makeResizable,removeDraggable,removeResizable,createMenuIcon,makeSortable,
    splitPane,getElDim,MDVProgress};