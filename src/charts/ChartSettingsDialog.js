import { BaseDialog } from "../utilities/Dialog.js";
import { createEl, createFilterElement } from "../utilities/Elements.js";
import noUiSlider from "nouislider";


/**
 * PJT: for review - this doesn't seem to be used anywhere?
 */
class ChartSettingsDialog extends BaseDialog{
    constructor(config,content){
        super(config,content);
        
    }

    init(content){
        
        for (let s of content){
            let d = createEl("div",{
                styles:{
                    padding:"5px"
                }
            },this.dialog);
            if (s.type !== "button"){
                createEl("label",{text:s.label},d);
            }
           
            if (s.info){
                this.addInfoIcon(s.info,d);
            }
            this[s.type](s,d);
        }

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
        const d1 =createEl("div",{},d)
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
            start: [s.current_value],
            range: {
                'min': [s.min],
                'max': [s.max]
            },
            step:s.step || null,
            tooltips:true,
            /*format:{
                to:v=>{
                    return v+""
                },
                from:v=>{
                    return Number(v)}
            },*/
            documentElement:s.doc
        });
        sl.noUiSlider.on("end",(values)=>{
            s.func(parseFloat(values[0]))
        })

    }

    check(s,d){
        const ch= createEl("input",{
            type:"checkbox",
        },d);
        ch.checked= s.current_value;
        ch.addEventListener("click",(e)=>{
            s.func(ch.checked)
        })
    }

    dropdown(s,d){
        const wrapper = createEl("div");
        const dd = createEl("select",{
            styles:{
                maxWidth:"200px"
            }
        }, wrapper);
        for (let item of s.values[0]){
            createEl("option",{
                text:item[s.values[1]],
                value:item[s.values[2]]
            },dd)
        }
        dd.value=s.current_value;
        dd.addEventListener("change",(e)=>{
            s.func(dd.value);
        });
        createFilterElement(dd, wrapper);
        d.append(wrapper);
    }

    doubleslider(s,d){
       
        const sl = createEl("div",{
            styles:{
              
            }
        },d)
        noUiSlider.create(sl, {
            start: [s.current_value[0], s.current_value[1]],
            range: {
                'min': [s.min],
                'max': [s.max]
            },
            documentElement:s.doc
        });
        sl.noUiSlider.on("end",(values)=>{
            s.func(parseFloat(values[0]),parseFloat(values[1]))
        })

    }

    text(s,d){
        const t = createEl("input",{
            value:s.current_value
        },d);
        t.addEventListener("keyup",()=>{
            s.func(t.value);
        })
    }

    textbox(s,d){
        const tb= createEl("textarea",{
           
            styles:{
                height:"20px",
                width:"100%"
            }
        },d);
        tb.value=s.current_value;
    
        tb.addEventListener("keypress",e=>{
            if (e.keyCode===13){
                tb.blur();
            }
            else{
                s.func(tb.value);
            }
        });
        tb.addEventListener("blur",e=>{
            tb.style.height = "20px";
            s.func(tb.value);
        })
        tb.addEventListener("focus",e=>{
            tb.style.height="80px"
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

export {ChartSettingsDialog}