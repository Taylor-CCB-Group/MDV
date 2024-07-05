onmessage=(e)=> {
    const {orderBuffer,filterBuffer,dimensions,data,displayData,colorScale,colorOnly}=e.data;
    //get access to the shared buffers
    const len = dimensions.length;
  
    const rowPositions = new Uint16Array(data,0,len);
    const colPositions = new Uint16Array(data,len*2,len);
    const values = new Uint8Array(data,len*4,len);

    const x= new Float32Array(displayData,0,len);
    const y= new Float32Array(displayData,len*4,len);
    const colors= new Uint8Array(displayData,len*8,len*3);
    const order= new Uint32Array(orderBuffer);
    const filter = new Uint8Array(filterBuffer);


    const c =  dimensions.columns;
    const currentRows = {};
    const reverse_map= {};
   
    let i=0;
   
    for (const n of order){
        if (filter[n]===0){
            currentRows[n]=i;
            reverse_map[i]=n;
            i++
        }
    }
    const cpg = c/dimensions.groups;

    i=0;

    const interval_size = (colorScale.max-colorScale.min)/colorScale.bins;
    
    for (let n=0;n<len;n++){
        const row = currentRows[rowPositions[n]]
        if (row===undefined){
            continue;
        }
        //update the colors
        let v = Math.min(colorScale.max,values[n]);
        v= Math.max(colorScale.min,v);   
        const color = colorScale.colors[Math.floor((v - colorScale.min) / interval_size)];    
        const p = 3*i;
        colors[p]=color[0];
        colors[p+1]=color[1];
        colors[p+2]=color[2];
        if (colorOnly){
            i++;
            continue;
        }
        const col = colPositions[n]; 
        x[i]=col*20 + (Math.floor(col/cpg)*200)+10;
        y[i]=row*20+10;
        i++;
    } 
    postMessage({reverse_map,total:i,currentRows});
    
}
