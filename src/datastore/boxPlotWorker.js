// biome-ignore lint/suspicious/noGlobalAssign: relatively innocuous in simple web worker
onmessage= function(e){
    if (e.data[4].analysis==="multi"){
        this.postMessage(multiBoxPlot(e));
        return;
    }
    const lFilter=new  Uint8Array(e.data[0]);
    const gFilter = new  Uint8Array(e.data[1]);
    const categories= new Uint8Array(e.data[2]);
    const arrType= e.data[3][1] === "int32"?Int32Array:Float32Array;
    const values = new arrType(e.data[3][0]);
    const config = e.data[4];
    const len =values.length;
    const catLen = config.values.length;

    const cats= new Array(catLen).fill(0);

    //get the number in each then create arrays
    //even though 2 iterations over data this is much faster
    for (let i=0;i<len;i++){
        if (gFilter[i]!==0){
            if  (gFilter[i] !==lFilter[i]){
            continue;
            }           
        }
        if (Number.isNaN(values[i])){
            continue;
        }        
        cats[categories[i]]++;
    }
    const arrs = [];
    for (let n=0;n<catLen;n++){
        arrs.push(new arrType(cats[n]))
    }

    const totals = new Array(catLen).fill(0);
    const xPos= new Float32Array(config.xPos);
 
    for (let n= 0;n<len;n++){
        if (gFilter[n]!==0){
            if  (gFilter[n] !==lFilter[n]){
            continue;
            }           
        }
        if (Number.isNaN(values[n])){
            continue;
        }  
        const cat = categories[n];
        arrs[cat][totals[cat]++]=values[n];
    }
    const boxStats=[];
 

    if (config.analysis==="boxplot"){
        const ids={};
        for (let n=0;n<catLen;n++){
            const arr= arrs[n];
            const l = arr.length
            if (l===0){
                continue;
            }
            arr.sort();
            const midIndex= median(0,l);
            const med = arr[midIndex]
            const Q1= arr[median(0,midIndex)];
            const Q3=arr[median(midIndex+1,l)];
            const IQR = Q3-Q1;
            let min=Q1-(1.5*IQR);
            let max=Q3+(1.5*IQR);
            min = min<arr[0]?arr[0]:min;
            max= max>arr[l-1]?arr[l-1]:max;
          
            ids[n]=boxStats.length; 
            boxStats.push({
                max:max,
                min:min,
                Q1:Q1,
                Q3:Q3,
                med:med,
                id:n
            });
            
        }
        for (let n= 0;n<len;n++){
            const p =  ids[categories[n]];  
            if (p!== undefined){
                xPos[n]= p*50+4+Math.random()*42;
            }
        }
    }
    else{
       
        const ranges = {};
        for (let n=0;n<catLen;n++){
            let k = null;
            if (config.scaletrim){
                k = stKernelDensityEstimator(arrs[n],config.ticks,config.scaletrim,config.bandwidth || 7);
            }
            else{
                k = kernelDensityEstimator(arrs[n],config.ticks,config.bandwidth || 7);
            }
            
            const max = Math.max(...k);
            
            if (max !==0 && !(Number.isNaN(max))){
                k.id=n;
                k.max=max;
                ranges[n]=[boxStats.length,max,k]
                boxStats.push(k);
               
            }
                     
        }
        const min = config.ticks[config.ticks.length-1];
        const interval =  config.ticks[0]-config.ticks[1];
        for (let n= 0;n<len;n++){
            
            const r = ranges[categories[n]];
            if (r){
                //const index =Math.floor((values[n]-min)/interval)
                //const fr = r[2][17-index]/r[1];
                xPos[n]= r[0]*50+4+Math.random()*42;
            }
        }
    }
    postMessage(boxStats);
}

function median(l,r){    
    return l+Math.round((r-l)/2)
}

function stKernelDensityEstimator(V,X,Qs,k=2){

    const den = Qs[1]-Qs[0];
    return X.map((x) => {
        let sum=0;
        let count=0;
        for (let i=0;i<V.length;i++){
            let val = V[i];
            if (!Number.isNaN(val)){
                val=(val-Qs[0])/den;
                val= val<0?0:val>1?1:val;
                let v= x-val;
                sum+=Math.abs(v /= k) <= 1 ? 0.75 * (1 - v * v) / k : 0;
                count++;
            }
        }
      return  sum/count;
    });
       
       
   
}

function kernelDensityEstimator(V,X,k=2) {   
    return X.map((x) => {
        let sum=0;
        let count=0;
        for (let i=0;i<V.length;i++){
            const val = V[i];
            if (!Number.isNaN(val)){
              let v= x-val;
              sum+=Math.abs(v /= k) <= 1 ? 0.75 * (1 - v * v) / k : 0;
              count++;
            }
        }
      return  sum/count;
    });
  
}


function multiBoxPlot(e){

    const data= [];
    for (const item of e.data[2]){
        data.push(item[1]==="int32"?new Int32Array(item[0]):new Float32Array(item[0]))
    }
    const dLen= data.length;
    const config = e.data[4];
   
    const scaleVals = config.scaleVals;
    for (const sv of scaleVals){
        sv.push(sv[1]-sv[0])
    }
  
    const lFilter=new  Uint8Array(e.data[0]);
    const gFilter = new  Uint8Array(e.data[1]);
    const tLen= gFilter.length;
    const results = new Array(dLen);
    const cats= new  Uint8Array(e.data[3]);
    const cat= config.cat;
    let total=0;
    for (let i=0;i<tLen;i++){
        if (cats[i] === cat && gFilter[i]===0){
             total++;  
        }         
    }

    for (let n=0;n<dLen;n++){
        //should cope with int32 and float32 
        results[n]=new Float64Array(total);
    }
    count=0;
    for (let i=0;i<tLen;i++){
        if (cats[i] === cat && gFilter[i]===0){     
            for (let n=0;n<dLen;n++){
                const v= (data[n][i]-scaleVals[n][0])/scaleVals[n][2];
                results[n][count]=v<0?0:v>1?1:v;   
            }
            count++;
        }
    }
    const boxStats=[];
    let tmin = Number.MAX_VALUE;
    let tmax = Number.MIN_VALUE
    for (let n=0;n<dLen;n++){
        const arr= results[n];
        const l = arr.length;  
        arr.sort();
        const midIndex= median(0,l);
        const med = arr[midIndex]
        const Q1= arr[median(0,midIndex)];
        const Q3=arr[median(midIndex+1,l)];
        const IQR = Q3-Q1;
        let min=Q1-(1.5*IQR);
        let max=Q3+(1.5*IQR);
        min = min<arr[0]?arr[0]:min;
        max= max>arr[l-1]?arr[l-1]:max;
        tmin= Math.min(tmin,min);
        tmax=Math.max(tmax,max);
        boxStats.push({
            max:max,
            min:min,
            Q1:Q1,
            Q3:Q3,
            med:med,
            id:n
        });   
    }
    boxStats.max=tmax;
    boxStats.min=tmin;
    return boxStats



}