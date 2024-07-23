
import BaseChart from "./BaseChart.js"
import ImagePanZoom from "../utilities/PanZoom.js";
import { createEl,makeSortable} from "../utilities/Elements.js";
import { getProjectURL } from "../dataloaders/DataLoaderUtil.ts";

class RowSummaryBox extends BaseChart{
    constructor(dataStore,div,config){
		super(dataStore,div,config);   
       
        if (config.image){
            const d = createEl("div",{
                classes:["mdv-image-holder"],
                style:{height:`${Number.parseInt(this.width*0.7)}px`}
            },this.contentDiv);       
            this.imViewer=new ImagePanZoom(d,this._getImage(0));
        }
        this.paramHolders={};
        const sectionDiv= createEl("div",{},this.contentDiv);
        for (const p of this.config.param){
            if (p === this.img_param){
                continue;
            }
            const c = this.dataStore.columnIndex[p];
            const d = createEl("div",{classes:["mdv-section"]},sectionDiv);
            d.__mcol__=c.field;
            createEl("div",{text:c.name,classes:["mdv-section-header"]},d);
            let te = createEl("div",{classes:["mdv-section-value"]},d);
            if (c.is_url){
                te= createEl("a",{target:"_blank"},te);
            }
            this.paramHolders[c.field]=te;
        }
        makeSortable(sectionDiv,{
            handle:"mdv-section-header",
            sortEnded:li=>{
                const c =  this.config
                let tp = [];
                if (c.image){
                    tp = [c.param[c.image.param]];
                    c.image.param=0;
                }
                c.param= tp.concat(li.map(x=>x.__mcol__));  
            }
        })
        this.setParam(0);
        this.contentDiv.style.overflowY="auto";
	}

    setParam(index){
       
        for (const p in this.paramHolders){
            const data  = this.dataStore.getRowText(index,p)
            this.paramHolders[p].textContent=data;
            if (this.paramHolders[p].nodeName==="A"){
                this.paramHolders[p].href=data;
            }
        }
    }

    _getImage(index){
        const c = this.config
        const  p  = c.param[c.image.param];
        this.img_param=p;
        const data =this.dataStore.getRowText(index,p);
        return getProjectURL(`${c.image.base_url}${data}.${c.image.type}`);

    }

    setSize(x,y){
        super.setSize(x,y);
        if (this.imViewer){
            this.imViewer.container.style.height= `${Math.round(this.width*0.7)}px`;
            this.imViewer.fit();
        }
    }

    onDataHighlighted(data){
        const i = data.indexes[0]
        if (this.imViewer){
            this.imViewer.setImage(this._getImage(i));
        }
        this.setParam(i);
    }

    onDataFiltered(){
         
    }

    changeBaseDocument(doc){
        super.changeBaseDocument(doc);
        this.imViewer.img.__doc__=doc;
    }

    getSettings(){
        return super.getSettings();     
    }
    
}

BaseChart.types["row_summary_box"]={
    "class":RowSummaryBox,
    name:"Row Summary Box",
    init:(config,dataSource,extraControls)=>{
        const is = extraControls.image_set
        if (is){
            if(is !== "__none__"){
                const li =dataSource.large_images[is];
                config.param = config.param || [];
                config.param.push(li.key_column);
                config.image={
                    base_url: getProjectURL(li.base_url),
                    type:li.type,
                    param:config.param.length-1
                }

            }  
        }
    },
    extra_controls:(dataSource)=>{
        //drop down of available image sets
        const li = dataSource.large_images;
        if (li){
            let vals=[];
            for (const n in li){
                vals.push({name:n,value:n});
            }
            vals = [{name:"none",value:"__none__"}].concat(vals);
            return [
                {
                    type:"dropdown",
                    name:"image_set",
                    label:"Image Set",
                    values:vals
                }
            ];
        }
        return []; 
    },
    params:[{
        type:"_multi_column:all",
        name:"Values To Display"
    }]
}

export default RowSummaryBox;