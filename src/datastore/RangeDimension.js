
import Dimension from "./Dimension.js";

class RangeDimension extends Dimension{
    /**
     * @param {DataStore} parent
     */
    constructor(parent) {
        super(parent);
        this.worker= new Worker (new URL("./binWorker.js?v=1",import.meta.url));  
    }

    filterSquare(args, columns) { 
        let data1= null;
        if (typeof columns[0] !== "string"){
            data1= columns[0];
        }
        else{
            data1 = this.parent.columnIndex[columns[0]].data;
        }
        
        const data2 = this.parent.columnIndex[columns[1]].data;
        const range1 = args.range1;
        const range2 = args.range2;
        const predicate = i => {
            const v1 = data1[i];
            const v2 = data2[i];
            return v1>=range1[0] && v1<=range1[1] && v2>=range2[0] && v2<=range2[1] && !isNaN(v1) && !isNaN(v2);
        }
        return this.filterPredicate({predicate}, columns);
    }

    /**
     * @param {Array<[number, number]>} args
     */
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
        const vs = points;

        const predicate = i => {
            const x = data1[i], y = data2[i];
            let inside = false;
            if (x<minX || x>maxX || y<minY || y>maxY || isNaN(x) || isNaN(y)){
                return false;
            }
            for (let i = 0, j = vs.length - 1; i < vs.length; j = i++) {
                let xi = vs[i][0], yi = vs[i][1];
                let xj = vs[j][0], yj = vs[j][1];

                let intersect = ((yi > y) != (yj > y))
                    && (x < (xj - xi) * (y - yi) / (yj - yi) + xi);
                if (intersect) inside = !inside;
            }
            return inside;
        }

        return this.filterPredicate({predicate}, columns);
    }

    filterRange(args,columns){
        const min = args.min;
        const max=args.max;   
        const arr = this.parent.columnIndex[columns[0]].data;
        // performance seems similar to non-predicate version
        const predicate = i => {
            const v = arr[i];
            return v >= min && v <= max && !isNaN(v);
        }
        return this.filterPredicate({predicate}, columns);
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
            [col.buffer,col.datatype],
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