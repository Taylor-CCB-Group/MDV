




onmessage= function(e){
    const config = e.data[3];
    const dtype = config.datatype=="text"?Uint8Array:Uint16Array
    const data=new dtype(e.data[2]);
    const lFilter=new  Uint8Array(e.data[0]);
    const gFilter = new  Uint8Array(e.data[1]);
    let result = null;
    if (config.method==="sankey"){
        const data2= new Uint8Array(e.data[4]);
        result = getSankeyData(lFilter,gFilter,data,data2,config)
    }
    else if (config.method==="venn"){
        result =calculateSets(lFilter,gFilter,data,config)
    }
    else if (config.method==="proportion"){
        const data2= new Uint8Array(e.data[4]);
        result = getProportionData(lFilter,gFilter,data,data2,config);
    }
    else if (config.method==="stacked"){
        const data2= new Uint8Array(e.data[4]);
        result = getStackedData(lFilter,gFilter,data,data2,config);
    }
    else if (config.method==="double_cat"){
        const data2 = config.datatype2==="multitext"?new Uint16Array(e.data(4)):new Uint8Array(e.data(4));
        result = getDoubleCategory(lFilter,gFilter,data,data2,config)
    }
    else{     
        result =getNumberInCategory(lFilter,gFilter,data,config)
    }
    postMessage(result);
}


function getDoubleCategory(lFilter,gFilter,data,data2,config){
    const len = data.length;
    const mt2 = config.datatype2==="multitext";
    const matrix = Array.from(config.values,(x,i1)=>{
        return Array.from(config.values2,(x,i2)=>({i1:i1,i2:i2,c:0}))       
    });
    for (let i=0;i<len;i++){
        if (gFilter[i]!==0){
            if  (gFilter[i] !==lFilter[i]){
                continue;
            }           
        }   
        if (mt2){
            const st = i*config.stringLength2;
            for (let n=st;n<config.stringLength2;n++){
                if (data[n]===65535){
                    break;
                }
                matrix[data[i]][data2[n]].c++;
            }
        }
        else{
            matrix[data[i]][data2[i]].c++
        }
    }
    return matrix;

}



//data the x category (groups) 
function getProportionData(lFilter,gFilter,data,data2,config){
    const len1 = config.values.length;
    const len2 = config.values2.length;
    const data3 = config.cats?new Uint8Array(config.cats):null;
    const len = data.length;
    const counts = new Array(len1);
    const totals = config.diviser?config.diviser:new Array(len1);
    for (let n=0;n<len1;n++){
        counts[n]=new Array(len2).fill(0);
        totals[n]=new Array(len2).fill(0);
    }
    const cat = config.category;
    let total=0;
    for (let i=0;i<len;i++){
        //if filtered out in global but not in local       
        totals[data[i]][data2[i]]++;    
        if ( gFilter[i]===0){
            if (data3 && data3[i]!==cat){
                continue;
            }
            counts[data[i]][data2[i]]++;
            total++;     
        }
    }
    let t_max=0;
    let t_min=10000000;
    for (let i =0;i<totals.length;i++){     
        const  t = totals[i];
        const c= counts[i];
        const nc=[];
        const vls = [];
        let total=0;
        let max=0;
        let min=10000000;
        for (let n=0;n<t.length;n++){
            if (t[n]===0){
                continue;
            }
            const v= config.denominators?c[n]/config.denominators[n]:(c[n]/t[n])*100;
            nc.push([v,i,n,Math.floor(Math.random()*6)]);
            vls.push(v);
            total+=v;
            max=Math.max(max,v);
            min = Math.min(min,v);
        }
        
        nc.av= total/nc.length;
        nc.std = std(vls,nc.av);
        nc.max=max;
        nc.min=min;
        counts[i]=nc
        t_max=Math.max(t_max,max);
        t_min = Math.min(t_min,min);
    }
    counts.max=t_max;
    counts.min=t_min;
    return counts;
}




function calculateSets(lFilter,gFilter,data,config){
    const t = performance.now();
    const len = data.length/config.stringLength;
    
    const sets = new Map();
    const ilen = config.stringLength;
    for (let i=0;i<len;i++){
        //if filtered out in global but not in local
        if (gFilter[i]!==0){
            if  (gFilter[i] !==lFilter[i]){
                continue;
            }           
        }
        const st = i*ilen;
        let arr = data.slice(st,st+ilen).toString();

        
        const val = sets.get(arr);
        if (!val){
            sets.set(arr,1);
        }
        else{
            sets.set(arr,val+1);
        }
        
        
    }
    const ret=[];
    for (const [key, value] of sets.entries(sets)) {
        const set_arr = key.split(",").map(x=>config.values[x]).filter(x=>x!==undefined);
        ret.push({sets:set_arr,size:value})
    }
    console.log(`calc ven : ${performance.now()-t}`);
    return ret;
   

    



}

function getNumberInCategory(lFilter,gFilter,data,config){
    
    const len = config.datatype==="multitext"?data.length/config.stringLength:data.length;
    const cats = new Array(config.values.length).fill(0);
    
    if (config.datatype==="multitext"){
        //calculateSets(lFilter,gFilter,data,config);
        const ilen = config.stringLength;
        for (let i=0;i<len;i++){
            //if filtered out in global but not in local
            if (gFilter[i]!==0){
                if  (gFilter[i] !==lFilter[i]){
                    continue;
                }           
            }
            const st = i*ilen;
            for (let n=st;n<st+ilen;n++){
                if (data[n]===65535){
                    break;
                }
                cats[data[n]]++;
            }
        }

    }
    else{
        for (let i=0;i<len;i++){
            //if filtered out in global but not in local
            if (gFilter[i]!==0){
                if  (gFilter[i] !==lFilter[i]){
                    continue;
                }           
            }
            cats[data[i]]++;
        }
    }
    return cats;

}

function getStackedData(lFilter,gFilter,data,data2,config){
    const len = data.length;
    const matrix = Array.from(config.values,(x,i)=>{
        return {
            id:i,
            values:Array.from(config.values2,(x,i)=>({id:i,count:0})),
            total:0
        }
    });
    for (let i=0;i<len;i++){
        if (gFilter[i]!==0){
            if  (gFilter[i] !==lFilter[i]){
                continue;
            }           
        }
        matrix[data[i]].values[data2[i]].count++;
        matrix[data[i]].total++;
    }
    for (let r of matrix){
        let rt = 0;
        let rpt=0;
        for (let i of r.values){
            const per= r.total===0?0:i.count/r.total;
            i.pos=rt;
            i.per=per;
            i.perpos=rpt;
            rt+=i.count;
            rpt+=per;
        }
    }
    return matrix;
}

function getSankeyData(lFilter,gFilter,data,data2,config){
    const len1 = config.values.length;
    const len2 = config.values2.length;
    const len = data.length;
    const matrix = new Array(len1);
    const nodes1= config.values.map((x,i)=>{
        return "A"+i;
    });
    const nodes2= config.values2.map((x,i)=>{
        return"B"+i;
    });



    for (let n=0;n<len1;n++){
        matrix[n]=new Array(len2).fill(0);
    }
    const links= [];
    let total=0;
    for (let i=0;i<len;i++){
        //if filtered out in global but not in local
        if (gFilter[i]!==0){
            if  (gFilter[i] !==lFilter[i]){
                continue;
            }           
        }
        matrix[data[i]][data2[i]]++;
        total++;
    }
    
    const nodes1Used=new Set();
    const nodes2Used=new Set()
    const minValue =Math.round(total/300);

    for (let i1 =0;i1<len1;i1++){
        for (let i2 = 0;i2<len2;i2++){
            const v = matrix[i1][i2];
            if (v!==0){
                nodes1Used.add(nodes1[i1]);
                nodes2Used.add(nodes2[i2]);
                links.push({source:nodes1[i1],target:nodes2[i2],value:v<minValue?minValue:v,trueValue:v})
            }
            
        }
    }
    const minNodes = Math.min(nodes1Used.size,nodes2Used.size);
    const n1 = Array.from(nodes1Used).map(x=>{
        return {id:x,ind:x.substring(1),param:0};
    })
    const n2 = Array.from(nodes2Used).map(x=>{
        return {id:x,ind:x.substring(1),param:1};
    })

    return {
        links:links,
        nodes:n1.concat(n2),
        minNodes:minNodes
    }
}

function std(arr,av){
    let n=arr.length-1;
    n=n===0?1:n;
    let std = arr.reduce((prev,cur)=>prev+(Math.pow(cur-av,2)),0);
    return Math.sqrt(std/n);

}