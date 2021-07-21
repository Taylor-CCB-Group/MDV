const func = function(e){
    
	const data= new Float32Array(e.data[2]);
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

self.onmessage= function(e){  
    self.postMessage(func(e));
}
export {func};

