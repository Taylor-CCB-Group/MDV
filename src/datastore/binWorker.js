//the following causes errors in jsdoc as the ] is nor parsed
//does not seem to affect the output
/**
 * @param {SharedArrayBuffer} e.data[0] - local filterBuffer
 * @param {SharedArrayBuffer} e.data[1] - global filterBuffer
 * @param {SharedArrayBuffer} e.data[2] - data
 * @param {Object} e.data[3]- config
 * @param {Number} e.data[3].bins - number of bins
 * @param {Number} e.data[3].min - min value
 * @param {Number} e.data[3].max - max value
 */
const func = (e)=> {
    
    const arrType = e.data[2][1]==="int32"?Int32Array:Float32Array;
	const data= new arrType(e.data[2][0]);
    const config = e.data[3];
    const interval_size = (config.max-config.min)/(config.bins);
    const max=config.max;
    const min =config.min;
    const len =data.length;
    const  histogram = new Array(config.bins+1).fill(0);
    const lFilter=new  Uint8Array(e.data[0]);
    const gFilter = new  Uint8Array(e.data[1]); 
    for (let i=0;i<len;i++){
        //if filtered out in global but not in local
        if (gFilter[i]!==0){
            if  (gFilter[i] !==lFilter[i]){
                continue;
            }           
        }
        let v= data[i];
      
        v=v>max?max:v<min?min:v;
        histogram[Math.floor((v - min) / interval_size)]++;
    }
    return histogram;
}

self.onmessage= (e)=> {  
    if (e.data.length === undefined || typeof e.data === "string") {
        // WordCloud internal setZeroTimeout calls window.postmessage...
        // seems ok to ignore like this, but not a very clean solution
        // ? there could be other instances of unexpected messages
        // console.warn("unexpected message to binWorker", e.data);
        return;
    }
    self.postMessage(func(e));
}
//export {func};

