import Dimension from "./Dimension.js";
class CatColDimension extends Dimension{
    constructor(parent){
        super(parent)
        this.worker= new Worker (new URL("./catColWorker.js?v=3",import.meta.url));      
    }

    getAverages(callback,columns,config={}){
        config.method=config.method || "mean";
        const t = performance.now();
        const cIndex= this.parent.columnIndex;
        this.worker.onmessage=(e)=>{
            console.log(`calc averages  ${performance.now()-t}`);
            callback(e.data);
        }

        config.values=cIndex[columns[0]].values;
        const colBuffers=[];
        for (let n=1;n<columns.length;n++){
           colBuffers.push(this.parent.columnIndex[columns[n]].buffer);
        }
        this.worker.postMessage([
            this.filterBuffer,
            this.parent.filterBuffer,
            cIndex[columns[0]].buffer,
            colBuffers,
            config
        ]);
     
       
    }

    destroy(){
        super.destroy();
        this.worker.terminate();
    }
}

Dimension.types["catcol_dimension"]=CatColDimension;

export default CatColDimension; 