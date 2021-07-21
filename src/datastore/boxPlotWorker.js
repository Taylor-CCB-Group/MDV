onmessage= function(e){
    const lFilter=new  Uint8Array(e.data[0]);
    const gFilter = new  Uint8Array(e.data[1]);
    const categories= new Uint8Array(e.data[2]);
    const values = new Float32Array(e.data[3]);
    const config = e.data[4];
    const len =values.length;
    const catLen = config.values.length;

    const cats= new Array(catLen).fill(0);

    for (let i=0;i<len;i++){
        if (gFilter[i]!==0){
            continue;
        }           
        cats[categories[i]]++;
    }
    const arrs = [];
    for (let n=0;n<catLen;n++){
        arrs.push(new Float32Array(cats[n]))
    }

    const totals = new Array(catLen).fill(0)
    for (let n= 0;n<len;n++){
        const cat = categories[n];
        arrs[cat][totals[cat]++]=values[n];
    }
    const boxStats=[]
    for (let n=0;n<catLen;n++){
        const arr= arrs[n];
        const l = arr.length
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
        boxStats.push({
            max:max,
            min:min,
            Q1:Q1,
            Q3:Q3,
            med:med
        })
        
    }
    postMessage(boxStats);
}

function median(l,r){    
    return l+Math.round((r-l)/2)
}