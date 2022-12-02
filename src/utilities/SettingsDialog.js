import { BaseDialog } from "./Dialog.js";
import { createEl } from "./Elements.js";
import noUiSlider from "nouislider";



class SettingsDialog extends BaseDialog{
    constructor(config,content){
        super(config,content);
        
    }

    init(content){
        this.controls=[]
        
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
            const factory = this[s.type];
            if (!factory) {
                console.warn(`SettingsDialog doesn't have a method '${s.type}'`);
                continue;
            }
            this.controls[s.label] = factory(s,d);
        }

    }
    /** TODO */
    folder(s, d) {
        const container = d;
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
        d.style.textAlign = "center";
        //maybe consider making this a 'button'
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
            //PJT::: this.config was undefined
            documentElement: s.doc
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


    multidropdown(s,d){
        const dd = createEl("select",{
            multiple:true,
            styles:{
                maxWidth:"200px",
                height:"100px"
            }
        });
        createEl("br",{},d);
        for (let item of s.values[0]){
            const v =item[s.values[2]];
            const args = {
                text:item[s.values[1]],
                value:v
            }
            const sel = s.current_value.indexOf(v) !== -1;
            if (sel){
                args.selected=true;
            }

            createEl("option",args,dd)
        }
        d.append(dd);
        createEl("br",{},d);
        const b = createEl("button",{
            classes:["ciview-button-sm"],
            text:"Change"
        },d)
        b.addEventListener("click",(e)=>{
            s.func(Array.from(dd.selectedOptions).map(x=>x.value));
        })
        return dd;
    }
       
    

    

    dropdown(s,d){
        const dd = createEl("select",{
            styles:{
                maxWidth:"200px"
            }
        });
        createEl("br",{},d);
        for (let item of s.values[0]){
            createEl("option",{
                text:item[s.values[1]],
                value:item[s.values[2]]
            },dd)
        }
        d.append(dd);
        dd.value=s.current_value;
        dd.addEventListener("change",(e)=>{
            s.func(dd.value);
            if (s.onchange){
                s.onchange(this.controls,dd.value);
            }
        });
        return dd;
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
        d.style.display = "flex";
        d.style.alignItems = "center";
        const t = createEl("input",{
            value:s.current_value,
            styles: {flex: 2}
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

export default SettingsDialog