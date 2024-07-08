import SortableDimension from "./SortableDimension.js";
import Dimension from "./Dimension.js";
export default class DeepToolsDimension extends SortableDimension {
    constructor(parent) {
        super(parent);
        this.deepToolsWorker= new Worker (new URL("./deepToolsWorker.js?v=1",import.meta.url)); 
    }


    setParameters({groups,columns,rows,length,displayData,data}){
        this.dimensions = {groups,columns,rows,length};
        this.displayData = displayData;
        this.data = data;
    }


    updateData(colorScale,colorOnly=false){
        return new Promise((resolve,reject)=>{
            this.deepToolsWorker.onmessage = (e) => {
                resolve(e.data);
            }
            this.deepToolsWorker.postMessage({
                displayData:this.displayData,
                data:this.data,
                filterBuffer:this.parent.filterBuffer,
                dimensions:this.dimensions,
                orderBuffer:this.orderBuffer,
                colorScale,
                colorOnly
            });
        });
    }

    destroy(){
        super.destroy();
        this.deepToolsWorker.terminate();
    }
}

Dimension.types["deeptools_dimension"]=DeepToolsDimension;