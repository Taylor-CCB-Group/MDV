import BaseChart from "./BaseChart.js";
import {createEl} from "../utilities/Elements.js";


class TextBoxChart extends BaseChart{
    constructor(dataStore,div,config){
		super(dataStore,div,config);       
        const c = this.config;
        c.text= c.text || "";
        this.para = createEl("div",{
            styles:{
                padding:"5px"
            }
        },this.contentDiv);
        //conditional import of function which uses the packages marked and sanitize
        //decreases size of entry module
        import('../utilities/MarkdownText').then(({ default:renderText })=>{
            this.render= renderText;
            this.para.innerHTML = this.render(c.text);
        })     
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
                this.para.innerHTML = this.render(x);
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
