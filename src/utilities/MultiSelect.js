import { createEl,addElProps,getElDim } from "./Elements";

class MultiSelectDropdown{
    constructor(div,options,config){
        this.config= {
            search:true,
            height:'15rem',
            maxShow:4,
            placeholder:'select',
            txtSelected:'selected',
            txtAll:'All',
            txtRemove: 'Remove',
            txtSearch:'search',
            doc:document,
            ...config
        }
        addElProps(div,{class:'multiselect-dropdown'})
        this.div=div;
        this.options=new Map();
      
        this.div.addEventListener('click',e=>this._handleClick(e));
        this._addDocListener();
      
        if (options){
            options.forEach(x=>this.addOption(x))
        }
        this.refresh();
    }

    changeDocument(newDoc){
        this._removeList();
        this._addDocListener();
        this.config.doc=newDoc;

    }

    _removeList(){
        if (this.listWrap){
            this.listWrap.remove();
            this.listWrap=null;
        }
    }

    _addDocListener(){
        const d =  this.config.doc;
        if (this.docListener){
           d.removeEventListener("click",this.docListener)
        }
        this.docListener= e=>{
            if (this.listWrap && !(this.listWrap.contains(e.target)) ) {
                this._removeList();
                this._callback();
            }
        }
        d.addEventListener("click",this.docListener);
    }

    getSelectedValues(){
        return Array.from(this.options.values()).filter(x=>x.selected).map(x=>x.value);
    }
    unSelectAll(){
        Array.from(this.options.values()).forEach(x=>x.selected=false);
        this.refresh();
    }

    setSelected(selected){
        selected.forEach(x=>{
            const o = this.options.get(x);
            o.selected= o.value === selected
        });
        this.refresh();
    }

    refresh(){
        this.div.querySelectorAll('span.optext, span.placeholder').forEach(t=>this.div.removeChild(t));
        const sels = Array.from(this.options.values()).filter(x=>x.selected);
        const nex= this.options.size-sels.length;
        if( nex< this.config.maxShow && nex !==0 && this.options.size>this.config.maxShow){
            const ex = Array.from(this.options.values()).filter(x=>!x.selected);
            const msg = ex.map(x=>x.text).join(",")
            this.div.appendChild(createEl('span',{class:['optext'],text:`${msg} excluded`}));  
        }
        else if(sels.length>(this.config.maxShow)){
            this.div.appendChild(createEl('span',{class:['optext'],text:`${sels.length} ${this.config.txtSelected}`}));          
          }
        else{
          sels.forEach(x=>{
            const c=createEl('span',{class:'optext',text:x.text});
            const rb= createEl('span',{class:'optdel',text:'🗙',title:this.config.txtRemove},c);
            rb.__rmvalue__=x.value;    
            this.div.appendChild(c);
          });
        }
        if(0===sels.length) this.div.appendChild(createEl('span',{class:'placeholder',text:this.config.placeholder}));
        this.div.querySelectorAll('span.optext, span.placeholder');
        if (this.listWrap){
            Array.from(this.list.childNodes).forEach(x=>{
                const v=this.options.get(x.__value__);
                if (v.selected){
                    x.classList.add("cheked")
                }
                else{
                    x.classList.remove("cheked")
                }
                x.childNodes[0].checked=v.selected;

            });
        }
    }
    remove(){
        this._removeList();
        this.config.doc.removeEventListener("click",this.docListener);
    }
    _handleSearch(){
        this.list.querySelectorAll(":scope div:not(.multiselect-dropdown-all-selector)").forEach(d=>{
            const txt=d.querySelector("label").innerText.toUpperCase();
            const v = this.search.value.toUpperCase();
            const does = txt.includes(v) || v ==="";
            if (does){
                console.log(txt);
            }
            d.style.display=does?'block':'none';
        });
    }
    _bulkChange(type){
        Array.from(this.options.values()).forEach(x=>x.selected=type==="All"?true:type==="None"?false:!x.selected);
        this.refresh();
    }
    _callback(){
        if (this.config.onchange){
            this.config.onchange(this.getSelectedValues())
        }  
    }
    _handleClick(e){
        const t =  e.target
        if (t.__rmvalue__){
            this.options.get(t.__rmvalue__).selected=false;
            this._callback();
            this.refresh();
        }
        else{
           this.createDropDown();
        }
        e.stopPropagation();
    }
    _handleListClick(e){
        const t =  e.target
      
        if (t.tagName==="BUTTON"){
            if (t.textContent==="OK"){
                this._removeList();
                this._callback();    
            }
            else{
                this._bulkChange(t.textContent);
            }
            e.stopPropagation();
        }
        else if (t.__value__){
            const opt = this.options.get(t.__value__);
            opt.selected=!opt.selected;
            this.refresh();
            e.stopPropagation();         
        }
    }

    createDropDown(){
        if (this.listWrap){
            return;
        }
     
        const dim = getElDim(this.div);
        this.listWrap =createEl('div',
            {
                classes:['multiselect-dropdown-list-wrapper'],
                style:{
                    top:`${dim.top+dim.height}px`,
                    left:`${dim.left}px`,
                    width:`${dim.width}px`,
                }
            },
            this.config.doc.body);
        this.list=createEl('div',{class:'multiselect-dropdown-list',style:{maxHeight:"300px"}},this.listWrap);   
        
        Array.from(this.options.values()).forEach(x=>{
            const op=createEl('div',{})
            const ic=createEl('input',{type:'checkbox'},op);
            ic.checked=x.selected;
            op.appendChild(ic);
            const l = createEl('label',{text:x.text},op);
            ic.__value__=l.__value__=op.__value__=x.value;
            this.list.appendChild(op);
        });

        const menu = createEl("div",{},this.listWrap)
        this.search = createEl('input',{class:['multiselect-dropdown-search'].concat([this.config.searchInput?.class??'form-control']),
                        style:{width:'95%',display:this.config.search?'block':'none'},
                        placeholder:this.config.txtSearch},menu);
        createEl("button",{text:"All",classes:["ciview-button"]},menu);
        createEl("button",{text:"None",classes:["ciview-button"]},menu);
        createEl("button",{text:"Toggle",classes:["ciview-button"]},menu);
        createEl("button",{text:"OK",classes:["ciview-button"]},menu);
        this.listWrap.addEventListener('click',e=>this._handleListClick(e));
        this.search.addEventListener("keyup",e=>this._handleSearch(e));
        //does list go off bottom of screen
        const ddim = getElDim(this.listWrap);
        const b_space= window.innerHeight-ddim.top//space at bottom
        const t_space = ddim.top;  //space at top
        if ( b_space>ddim.height){
            //do nothing
        }else if (t_space>ddim.height){
            this.listWrap.style.top=`${dim.top-ddim.height}px`;
        }
        else{
            const m = Math.max(b_space,t_space)-60;
            this.list.style.height= `${m}px`;
            if (t_space>b_space){
                this.listWrap.style.top="10px";
            }
        }
       
      
    }
    addOption({text,value,selected}){
        this.options.set(value,{
            text:text,
            selected:selected,
            value:value
        });  
    }

    removeOption(value){
        this.options.delete(value);
        const el =Array.from(this.list.childNodes).find(x=>x.__value__===value);
        if (el){
            el.remove();
            return true;
        }
        return false;
    }
}

export default MultiSelectDropdown;