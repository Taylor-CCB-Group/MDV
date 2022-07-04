import {DataModel} from "../table/DataModel.js";
import BaseChart from "./BaseChart.js"
import { createEl } from "../utilities/Elements.js";
import ImageTable from "../table/ImageTable.js";



class ImageTableChart extends BaseChart{
    constructor(dataStore,div,config){
		super(dataStore,div,config);
       
        this.dataModel= new DataModel(dataStore,{
            autoupdate:false
        });
        this.dataModel.setColumns(this.config.param);
        this.dataModel.updateModel();
        const c = this.config;

     

        this.grid= new ImageTable(this.contentDiv,this.dataModel,{
            base_url:c.images.base_url,
            image_type:c.images.type,
            image_key:c.param[0],
            initial_image_width:c.image_width
        });
        this.grid.addListener("image_clicked",(e,index)=>{
            this.dataStore.dataHighlighted([index],this)
        },c.id);
      
        
	}

    setSize(x,y){
        super.setSize(x,y);
        this.grid.resize();
    }

    onDataFiltered(){
        this.dataModel.updateModel();
        this.grid.show();       
    }

    getColorOptions(){
        return {
            colorby:"all"
        }
    }

    onDataHighlighted(data){
        if (data.source===this){
            return;
        }
        const id = data.indexes[0];
        this.grid.scrollToTile(id,true);

    }

    colorByColumn(column){
		this.grid.setColorBy(this.getColorFunction(column,false));
        const ft = this.grid.getFirstTileInView();
        this.grid.show(ft);
	}

    getSettings(){
        const od = this.grid.originalDimensions;
        const c = this.config;
        const settings = super.getSettings();
        return settings.concat([{
            type:"slider",
            max:od[1]*2,
            min:10,
            doc:this.__doc__,
            label:"Image Size",
            current_value:c.image_width || od[0],
            func:x=>{
                c.image_width=x;
                this.grid.setImageWidth(x,true)
            }
        }])
    }
    changeBaseDocument(doc){
        super.changeBaseDocument(doc);
        this.grid.__doc__=doc;
      }
}



BaseChart.types["image_table_chart"]={
    "class":ImageTableChart,
    name:"Image Table",
    required:"images",
    init:(config,dataSource)=>{
        //set the chr,start,finish columns
        const i = dataSource.images;
        config.param= i.key
        //set the base track corresponding to the datastore
        config.images={
            base_url:i.base_url,
            type:i.type
        }
    }

}

export default ImageTableChart;