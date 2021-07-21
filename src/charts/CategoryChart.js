import {SVGChart} from "./SVGChart.js";

class CategoryChart extends SVGChart{
    constructor(dataStore,div,config,axisTypes){
        
        if (Array.isArray(config.param)){
            config.param=config.param[0];
        }
        config.title = config.title || dataStore.getColumnName(config.param);        
		super(dataStore,div,config,axisTypes);        		 
        this.dim = this.dataStore.getDimension("category_dimension",config.param);
        this.colors  = this.dataStore.getColumnColors(config.param);  
        this.filter=[];
       
     
	}

    remove(){
        this.dim.destroy();
        super.remove();
    }

    setSize(x,y){
        super.setSize(x,y);
        this.drawChart();
    }

    removeFilter(){
        this.dim.removeFilter();
        this.filter=[];  
        this.drawChart();   
    }

    filterCategories(cat,append){
        if (append){
            this.filter.push(cat)
        }
        else{
            this.filter=[cat];
        }
        this.resetButton.style.display = "inline";
        this.drawChart(100);
        this.dim.filterCategories(this.filter.lenght===1?this.filter[0]:this.filter);       
    }

    getFilter(){
        return this.filter.slice(0)
    }

    pinChart(){
        this.isPinned=true;
    }

    unpinChart(){
        this.isPinned = false;
        this.onDataFiltered();
    }

    onDataFiltered(dim){
        //no need to change anything
        if (this.dim === dim || this.isPinned){
            return;
        }
        if (dim === "all_removed"){
            this.filter=[];
            this.resetButton.style.display="none";
        }
        const config={};
      
        this.dim.getCategories(data=>{
            this.rowData=[];
            this.maxCount=1;
            for (let n=0;n<data.length;n++){
                this.rowData.push([data[n],n]);
                this.maxCount= Math.max(data[n],this.maxCount);
            }
            this.drawChart();            
        },config)    
    }
}

export {CategoryChart};
