import { BaseDialog } from "./Dialog.js";
import { createEl, createFilterElement } from "./Elements.js";
import noUiSlider, { create } from "nouislider";
// import lgui from 'lil-gui';
import { action } from "mobx";


class SettingsDialog extends BaseDialog{
    constructor(config,content){
        super(config,content);
    }
    
    init(content, parent){
        if (!this.controls) this.controls=[];
        if (!parent) parent = this.dialog;
        //experimental lil-gui version...
        // this.initLilGui(content);
        // return;
        for (let s of content){
            let d = createEl("div",{
                styles:{
                    // padding:"5px"
                }
            },parent);
            if (s.type !== "button"){
                createEl("label",{text:s.label},d);
            }
           
            if (s.info){
                this.addInfoIcon(s.info,d);
            }
            //not sure  why assign a varaiable
            //'this' no longer available in function call - causing problems
            //const factory = this[s.type];
            if (!this[s.type]) {
                console.warn(`SettingsDialog doesn't have a method '${s.type}'`);
                continue;
            }
            if (this.config.useMobx && s.func) { //mobx complains about mutating state outside of actions
                // we could check whether we'd already wrapped the function, but getSettings() returns a new object each time
                s.func = action(`${s.label} <action>`, s.func);
            }
            this.controls[s.label] = this[s.type](s,d); //this is going to go wrong if the label isn't unique - like if we have multiple similar layers
        }
    }
    initLilGui(content){
        const gui = new lgui({container: this.dialog});
        //gui.name = this.config.title + " Settings";
        for (let s of content) {
            switch (s.type) {
                case "dropdown":
                    gui.add(s, 'current_value', s.values[0].map(c => c[s.values[1]])).name(s.label).onChange(s.func);
                    break;
                case "check":
                    if (s.current_value === undefined) s.current_value = false;
                    gui.add(s, 'current_value').name(s.label).onChange(s.func);
                    break;
                case "button":
                    gui.add(s, 'func').name(s.label);
                    break;
                case "slider":
                    gui.add(s, 'current_value', s.min, s.max).name(s.label).onChange(s.func);
                    break;
                case "spinner":
                    gui.add(s, 'current_value', s.min, s.max).name(s.label).onChange(s.func);
                    break;
                case "radiobuttons":
                    gui.add(s, 'current_value', s.choices).name(s.label).onChange(s.func);
                    break;
                case "text":
                case "textbox":
                default:
                    gui.add(s, 'current_value').name(s.label).onChange(s.func);
            }
        }
    }
            
    /** TODO */
    folder(s, d) {
        this.init(s.current_value);
    }
    spinner(s,d){
         
        let sp =createEl("input",{
            type:"number",
            value:s.current_value,
            max:s.max || null,
            min:s.min || 0,
            step:s.step || 1
        },d);
        sp.addEventListener("change",(e)=>{
            s.func(parseInt(sp.value));
        });
    }
    radiobuttons(s,d){
        // createEl("br",{},d);
        const d1 =createEl("div",{
            classes:["ciview-radio-group"]
        },d)
        for (let c of s.choices){
            createEl("span",{   
                text:c[0]
            },d1);
            const rb = createEl("input",{   
                type:"radio",
                name:s.label,
                value:c[1]
            },d1);
            rb.checked = s.current_value===c[1];
            rb.addEventListener("click",(e)=>{
                s.func(e.currentTarget.value);
            })

        }
    }

    button(s,d){
        d.style.textAlign = "center";
        createEl("button",{
            classes:["ciview-button"],
            text:s.label
        },d)
        .addEventListener("click",(e)=>{
            s.func()
        })
    }

    slider(s,d){
        const sl = createEl("div",{},d)
        noUiSlider.create(sl, {
            start: [s.current_value || 0],
            range: {
                'min': [s.min || 0],
                'max': [s.max || 1]
            },
            step:s.step || null,
            tooltips:true,
            //PJT::: this.config was undefined (no longer the case)
            documentElement: s.doc || this.config.doc
        });
        const change = (values) => {
            s.func(parseFloat(values[0]))
        };
        sl.noUiSlider.on("end", change);
        if (s.continuous) {
            console.log('adding update listener for ', s.label);
            sl.noUiSlider.on("update", change);
        }
        return sl;
    }

    check(s,d){
        const ch= createEl("input",{
            type:"checkbox"
        },d);
        ch.checked= s.current_value;
        ch.addEventListener("click",(e)=>{
            s.func(ch.checked)
        })
    }


    multidropdown(s, d){
        const wrapper = createEl("div");
        const dd = createEl("select",{
            multiple:true,
            styles:{
                maxWidth:"200px",
                height:"100px"
            }
        }, wrapper);
        createEl("br",{},d);
        for (let item of s.values[0]){
            // pass 1d array for simple list of strings, or [object[], text_field, value_field] for objects
            const text = s.values.length > 1 ? item[s.values[1]] : item;
            const value = s.values.length > 1 ? item[s.values[2]] : item;
            const args = {
                text,
                value
            }
            const sel = s.current_value.indexOf(value) !== -1;
            if (sel){
                args.selected=true;
            }

            createEl("option",args,dd)
        }
        createEl("br",{},d);
        const b = createEl("button",{
            classes:["ciview-button-sm"],
            text:"Change"
        },d)
        b.addEventListener("click",(e)=>{
            s.func(Array.from(dd.selectedOptions).map(x=>x.value));
        })
        createFilterElement(dd, wrapper);
        d.append(wrapper);
        return dd;
    }
       
    

    

    dropdown(s, d){
        const wrapper = createEl("div");
        const dd = createEl("select",{
            styles:{
                maxWidth:"200px"
            }
        }, wrapper);
        // createEl("br",{},d);
        for (let item of s.values[0]) {
            // pass 1d array for simple list of strings, or [object[], text_field, value_field] for objects
            const text = s.values.length > 1 ? item[s.values[1]] : item;
            const value = s.values.length > 1 ? item[s.values[2]] : item;
            createEl("option",{
                text,
                value
            },dd)
        }
        dd.value=s.current_value;
        dd.addEventListener("change",(e)=>{
            s.func(dd.value,this.controls);
            dd.title = dd.value; // for tooltip, not sure if best accessibility practice
            if (s.onchange){
                s.onchange(this.controls,dd.value);
            }
        });
        createFilterElement(dd, wrapper);
        d.append(wrapper);
        return wrapper;
    }

    doubleslider(s,d){
       
        const sl = createEl("div",{
            styles:{
              
            }
        },d)
        noUiSlider.create(sl, {
            tooltips:true,
            start: [s.current_value[0], s.current_value[1]],
            range: {
                'min': [s.min],
                'max': [s.max]
            },
            documentElement:s.doc
        });
        const change = (values) => {
            s.func(parseFloat(values[0]), parseFloat(values[1]))
        };
        sl.noUiSlider.on("end", change);
        if (s.continuous) {
            sl.noUiSlider.on("update", change);
        }
        return sl;
    }

    text(s,d){
        // d.style.display = "flex";
        // d.style.alignItems = "center";
        const t = createEl("input",{
            value:s.current_value,
            size: 10,
            styles: {flex: 1}
        },d);
        t.addEventListener("keyup",(e)=>{
            if (s.only_update_on_enter){
                if (e.key==="Enter"){
                    s.func(t.value);
                }
                return;
            }
            s.func(t.value);
        })
    }

    textbox(s,d){
        const tb= createEl("textarea",{
           
            // no other 'textbox' in our code apart from TextBoxChart AFAICT.
            styles:{
                // height:"20px",
                width:"100%"
            }
        },d);
        tb.value=s.current_value;
    
        tb.addEventListener("input",e=>{
            s.func(tb.value);
        });
        tb.addEventListener("blur",e=>{
            s.func(tb.value);
        })
    }

    addInfoIcon(info,parent){
        
        const sp= createEl("span",{
            "aria-label":info.text,
            role:"tooltip",
            "data-microtip-size": "medium",
            "data-microtip-position": info.position,
            styles:{
                marginLeft:"5px"
            },
            classes:["ttt"]
            },parent);
    
        createEl("i",{  
            classes:["fas","fa-info"]
        },sp);
        return sp;
    }

}

export default SettingsDialog