import {BaseDialog} from "../../utilities/Dialog.js";
import {createEl} from "../../utilities/Elements.js";
import noUiSlider from "nouislider";


class CustomDialog extends BaseDialog{
    constructor(config){
        if (config.buttons){
            config.footer=true;
            config._buttons= config.buttons;
            config.buttons = undefined
        }
        config.width= config.width || 300;
        
        super(config,config);

    }

    init(config){
        this.controlValues={};
        this.controls={};
        if (config._buttons){
            for (const b of config._buttons){
                this.addButton(b);
            }
        }
        if (config.text){
            const t = createEl("div",{},this.dialog);
            t.innerHTML=config.text;
        }
     
        if (config.controls){
            for (const s of config.controls){
                if (!s.id){
                    s.id=s.label;
                }
                const d = createEl("div",{
                    styles:{
                        padding:"5px"
                    }
                },this.dialog);
                if (s.type !== "button"){
                    const dis = s.type==="checkbox"?"inline":"block"
                    createEl("label",{text:s.label,styles:{fontWeight:"bold",display:dis,marginBottom:"3px"}},d);
                }             
                this[s.type](s,d);               
                if (s.description){
                   
                    createEl("div",{text:s.description,styles:{fontSize:"12px"}},d);
                }
            }
        }
    }

    radiobuttons(s,d){
        const d1 =createEl("div",{},d);
        this.controlValues[s.id]=s.current_value;
        for (const c of s.choices){
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
                const v =  e.currentTarget.value;
                this.controlValues[s.id]=v;
                if (s.func){
                    s.func(v,this);
                }            
            })

        }
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
            documentElement:s.doc || document
        });
        this.controlValues[s.id]=s.current_value;
        sl.noUiSlider.on("end",(values)=>{
            const v= Number.parseFloat(values[0]);
            this.controlValues[s.id]=v;
            if (s.func){
                s.func(v);
            }
        });
    }

    changeDropDownContent(dropdown_id,items){
        const dd = this.controls[dropdown_id];
        while(dd.length!==0){
            dd.remove(0);
        }
        for (const item of items){
            createEl("option",{
                text:item.text,
                value:item.value
            },dd)
        }
        this.controlValues[dropdown_id]=dd.value;
    }


    dropdown(s,d){
        const dd = createEl("select",{
            styles:{
                maxWidth:"200px"
            }
        });
        for (const item of s.items){
            createEl("option",{
                text:item.text,
                value:item.value
            },dd)
        }
        d.append(dd);
        if (s.current_value){
            dd.value=s.current_value;
            
        }
        this.controlValues[s.id]=dd.value;
        dd.addEventListener("change",(e)=>{
            this.controlValues[s.id]=dd.value;
            if (s.func){
                s.func(dd.value,this);
            }
        });
        this.controls[s.id]=dd;
    }

    checkboxgroup(s,d){
        const self =this;
        const menu = createEl("div",{},d);
        function addButton(title,check){
            const sall = createEl("button",{
                text:title,
                classes:["ciview-button-sm"]
            },menu);
            sall.addEventListener("click",()=>{
                for (const ch of s.checkboxes){
                    self.controls[ch.id].checked=check;
                    self.controlValues[ch.id]=check;
                }
            });
        }

        addButton("Check All",true);
        addButton("Uncheck All",false);
       
       
        const cbHolder = createEl("div",{},d);
        for (const cb of s.checkboxes){
            const cbd =  createEl("div",{styles:{
                display:"inline-block",
                padding:"3px"
            }},cbHolder)
            createEl("span",{
                styles:{
                    fontSize:"13px"
                },
                text:cb.label

            },cbd);
            cb.id =`${s.id}_${cb.id}`
            this.checkbox(cb,cbd)
        }

    }

    checkbox(s,d){
        const sv = s.current_value;
        const ch= createEl("input",{
            type:"checkbox",
            classes:["mdv-checkbox"]
        },d);
        ch.checked=s.current_value;
        this.controlValues[s.id]=sv;
        ch.addEventListener("click",(e)=>{
            this.controlValues[s.id]=ch.checked;
            if (s.func){
                s.func(ch.checked);
            }
        });
        this.controls[s.id]=ch;
    }

    text(s,d){
     
        const v= s.defaultValue || "";
        const t = createEl("input",{
            value:v
        },d);
    
        this.controlValues[s.id]=v;
        t.addEventListener("keyup",()=>{
            this.controlValues[s.id]=t.value;
        });
        this.controls[s.id]=t;
    }


    

    buttonClicked(button){
        const rv = button.method(this.controlValues, this.controls);
        if (rv !== "noclose"){
            this.close();
        }
       
    }

    addButton(button){
        createEl("button",{
            text:button.text,
            classes:["ciview-button"]
        },this.footer).addEventListener("click",()=>this.buttonClicked(button));

    }
}


export default CustomDialog;