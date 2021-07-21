
onmessage= function(e){
    //0 means it is included
	const data= new Uint8Array(e.data[2]);
    const config = e.data[3];
    const  len = data.length;
    const cats = new Array(config.values.length).fill(0);
    const lFilter=new  Uint8Array(e.data[0]);
    const gFilter = new  Uint8Array(e.data[1]);
    for (let i=0;i<len;i++){
        //if filtered out in global but not in local
        if (gFilter[i]!==0){
            if  (gFilter[i] !==lFilter[i]){
            continue;
            }           
        }
        cats[data[i]]++;
    }
    postMessage(cats);
}