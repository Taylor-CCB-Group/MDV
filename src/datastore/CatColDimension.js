import Dimension from "./Dimension.js";
class CatColDimension extends Dimension{
    constructor(parent){
        super(parent)
        this.worker= new Worker (new URL("./catColWorker.js?v=7",import.meta.url));      
    }


    filterCatCol(args,columns){
        const p = this.parent;
        const filter = p.filterArray;
        const localFilter= this.filterArray;
        const catCol= p.columnIndex[columns[0]]
        const catData = catCol.data;
        const len = catData.length;
        const colData = p.columnIndex[columns[1]].data;
        const thr = args.threshold || 0;
        const cat = catCol.values.indexOf(args.cat);
     
        for (let i=0;i<len;i++){
            //repetitive code - ideally have inline function
            const d= colData[i];
            if (catData[i]===cat && !(Number.isNaN(d)) && d>thr){
                //is filtered locally
                if (localFilter[i]===1){
                    //remove from global filter
                    if(--filter[i]===0){
                        //increase unfiltered size if no filters left on this row
                        p.filterSize++;                    
                    }                   
                }
                //remove local filter
                localFilter[i]=0;
            }
            else{
                 //not already filtered locally, update global
                 if (localFilter[i]===0){
                    if(++filter[i]===1){
                        p.filterSize--;                   
                    }
                }
                //add local filter
                localFilter[i]=1;
            }
        }
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
        config.columns=columns.slice(1);
        const colBuffers=[];
        for (let n=1;n<columns.length;n++){
           colBuffers.push([cIndex[columns[n]].buffer,cIndex[columns[n]].datatype]);
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