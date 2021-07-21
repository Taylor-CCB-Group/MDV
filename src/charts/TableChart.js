import { SlickGrid } from "../table/SlickGrid.js";
import { RowSelectionModel } from "../table/RowSelectionModel.js";
import {DataModel} from "../table/DataModel.js";
import {BaseChart} from "./BaseChart.js"


class TableChart extends BaseChart{
    constructor(dataStore,div,config){
		super(dataStore,div,config);  
        let cols = dataStore.getColumnList();
        cols = [{field:"index",id:"index",name:"index"}].concat(cols);
        for (let col of cols){
            col.id= col.field;
            col.sortable=true;
        }
       this.options = {
            enableCellNavigation: true,
            enableColumnReorder: false
        };	
        this.dataModel= new DataModel(dataStore,{autoupdate:false});
        this.grid= new SlickGrid(this.contentDiv,this.dataModel,cols,this.options);
        this.dataModel.addListener(this.config.id,()=>{
            this.grid.invalidate();
        });
        this.grid.setSelectionModel(new RowSelectionModel());
        this.grid.init();


        this.grid.onSort.subscribe( (e, args)=> {
            this.dataModel.sort(args.columnId,args.sortAsc?"asc":"desc");
            this.grid.invalidateAllRows();
            this.grid.render();
        });
        this.filter=[];
        this.onDataFiltered();
     
	}


    changeBaseDocument(doc){
        super.changeBaseDocument(doc);
        this.grid.setBaseDocument(doc);
    
    }

    setSize(x,y){
        super.setSize(x,y);
        this.grid.resizeCanvas();
        this.grid.invalidate();
    }
 
    pinChart(){
        this.isPinned=true;
    }

    unpinChart(){
        this.isPinned = false;
        this.onDataFiltered();
    }

    addColumns(){

    }

    onDataFiltered(){
        if (this.isPinned){
            return;
        }
        this.dataModel.updateModel()
    }

    
}

export {TableChart}

BaseChart.types["table_chart"]={
    "class":TableChart,
    name:"Table",
    params:[]

}