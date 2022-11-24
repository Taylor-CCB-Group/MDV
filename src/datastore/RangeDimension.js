
import Dimension from "./Dimension.js";

class RangeDimension extends Dimension{
    constructor(column,parent){
        super(column, parent);   
        this.worker= new Worker (new URL("./binWorker.js",import.meta.url));  
    }

    filterSquare(args,columns){
     
        let data1= null;
        if (typeof columns[0] !== "string"){
            data1= columns[0];
        }
        else{
            data1 = this.parent.columnIndex[columns[0]].data;
        }
        
        const data2 = this.parent.columnIndex[columns[1]].data;
        const filter = this.parent.filterArray;
        const parent = this.parent;
        const range1 =args.range1;
        const range2 =args.range2;
      
        const localFilter= this.filterArray;
        for (let i=0;i<this.parent.size;i++){
            const v1 = data1[i];
            const v2 = data2[i];
            if ( v1<range1[0] || v1>range1[1] || v2<range2[0] || v2>range2[1] || isNaN(v1) || isNaN(v2)){
               
                if (localFilter[i]===0){
                    if(++filter[i]===1){
                        parent.filterSize--;
                                         
                    }
                }                  
                localFilter[i]=1
            }
            else{
                if (localFilter[i]===1){
                    if(--filter[i]===0){
                        parent.filterSize++;                    
                    }                   
                }
                localFilter[i]=0;
            }
        }
    }

    filterPoly(args,columns){
        const points=args;
        let minX=Number.MAX_VALUE, minY= Number.MAX_VALUE;
        let maxX=Number.MIN_VALUE, maxY= Number.MIN_VALUE;
        for (let pt of points){
            minX=Math.min(minX,pt[0]);
            maxX= Math.max(maxX,pt[0]);
            minY=Math.min(minY,pt[1]);
            maxY= Math.max(maxY,pt[1]);

        }
        let data1= null;
        if (typeof columns[0] !== "string"){
            data1= columns[0];
        }
        else{
            data1 = this.parent.columnIndex[columns[0]].data;
        }
        const data2 = this.parent.columnIndex[columns[1]].data;
        const filter = this.parent.filterArray;
        const parent = this.parent;
        const len = parent.size;
        const vs =points;
        const localFilter= this.filterArray;
        for (let n=0;n<len;n++){
			let x = data1[n], y = data2[n];
			let inside = false;
            if (x<minX || x>maxX || y<minY || y>maxY || isNaN(x) || isNaN(y)){
                if (localFilter[n]===0){
                    if(++filter[n]===1){
                        parent.filterSize--;
                                         
                    };
                }                  
                localFilter[n]=1
            }
            else{
                for (let i = 0, j = vs.length - 1; i < vs.length; j = i++) {
                    let xi = vs[i][0], yi = vs[i][1];
                    let xj = vs[j][0], yj = vs[j][1];

                    let intersect = ((yi > y) != (yj > y))
                        && (x < (xj - xi) * (y - yi) / (yj - yi) + xi);
                    if (intersect) inside = !inside;
                }

                if(!inside){
                    if (localFilter[n]===0){
                        if(++filter[n]===1){
                            parent.filterSize--;
                                            
                        };
                    }                  
                    localFilter[n]=1;
                }
                else{
                    if (localFilter[n]===1){
                        if(--filter[n]===0){
                            parent.filterSize++;                    
                        }                   
                    }
                    localFilter[n]=0;
                }
            }		
		}
    }

    filterRange(args,columns){
        const min = args.min;
        const max=args.max;   
        const arr = this.parent.columnIndex[columns[0]].data;
        const filter = this.parent.filterArray;
        const localFilter= this.filterArray;
        const parent = this.parent;
        for (let i=0;i<this.parent.size;i++){
            const v = arr[i];
            if (v<min || v>max || isNaN(v)){
                if (localFilter[i]===0){
                    if(++filter[i]===1){
                        parent.filterSize--;
                    };
                }             
                localFilter[i]=1
            }
            else{
                if (localFilter[i]===1){
                    if(--filter[i]===0){
                        parent.filterSize++;
                    }
                }
                localFilter[i]=0;
            }
        }   
    }

    getBins(callback,column,config={}){
       
        const col = this.parent.columnIndex[column];
        config.bins = config.bins===undefined?10:config.bins;
        config.min = config.min===undefined?col.minMax[0]:config.min;
        config.max= config.max===undefined?col.minMax[1]:config.max;
       
        
        const t = performance.now();
        const action  =(e)=>{
            console.log(`calculate bins ${col.name} : ${performance.now()-t}`);
            callback(e.data);
                
        };
        this.worker.onmessage=action;
           
        this.worker.postMessage([
            this.filterBuffer,
            this.parent.filterBuffer,
            col.buffer,
            config
            ]);
    }

    destroy(){
        super.destroy();
        this.worker.terminate();
    }
}

Dimension.types["range_dimension"]=RangeDimension;

export default RangeDimension;