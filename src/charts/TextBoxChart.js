import BaseChart from "./BaseChart.js";
import {createEl} from "../utilities/Elements.js";
import {marked} from "marked";
import {sanitize} from 'dompurify';

function render(text) {
    return sanitize(marked.parse(text));
}

class TextBoxChart extends BaseChart{
    constructor(dataStore,div,config){
		super(dataStore,div,config);
       
        const c = this.config;
        c.text= c.text || "";
        this.para = createEl("div",{
            // text: c.text,
            // innerHTML: marked.parse(c.text),
            styles:{
                padding:"5px"
            }
        },this.contentDiv);
        this.para.innerHTML = render(c.text);

        
	}

    getSettings(){
        const settings = super.getSettings();
        const c= this.config;
        settings.push({
            label:"Text",
            type:"textbox",
            current_value:c.text,
            func: x => {
                c.text = x;
                this.para.innerHTML = render(x);
                //this.para.textContent = c.text = x
            }
        });
        return settings;
    }

}

export default TextBoxChart;

BaseChart.types["text_box_chart"]={
    "class":TextBoxChart,
    name:"Text Box",
    params:[],
    extra_controls: ()=> [
        {
            type:"textbox",
            name:"text",
        }
    ],
    init:(config, dataSource, extraControls) => {
        config.text = extraControls.text;
    }
}
