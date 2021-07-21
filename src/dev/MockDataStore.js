import DataStore from "../datastore/DataStore.js";

function getRandomDataStore(size,config={}){
    const columns= [
        {
            name:"x1",
            field:"x1",
            datatype:"double",
        },
        {
            name:"x2",
            field:"x2",
            datatype:"double"
        },
        {
            name:"x3",
            field:"x3",
            datatype:"double"
        },
        {
            name:"colors",
            field:"x4",
            datatype:"text",
            values:[
            "blue", 
            "purple",
            "red",
            "yelllow"
            ],
            colors:[
            "#7FFFD4",
            "#8A2BE2",
            "#DC143C",
            "#FFD700"
            ]
        }
        
    ];
    const c1 = new SharedArrayBuffer(size*4);
    const c2 = new SharedArrayBuffer(size*4);
    const c3 = new SharedArrayBuffer(size*4);
    const c4 = new SharedArrayBuffer(size*1);

    const a1 = new Float32Array(c1);
    const a2 = new  Float32Array(c2);
    const a3 = new  Float32Array(c3);
    const a4 = new Uint8Array(c4);


    config.ranges=config.ranges || [[-50,50],[-50,50],[-50,50]]
    const cols = [a1,a2,a3];
    const m_factors =config.ranges.map(x=>[x[1]-x[0],x[0]])


    for (let n=0;n<size;n++){
        for (let i=0;i<3;i++){
            cols[i][n]=(Math.random()*m_factors[i][0])+m_factors[i][1];
        }
        
        a4[n]=Math.floor((Math.random()*4));
    
        if (Math.random()>0.95){
            a1[n]=NaN;
        }
    
    }
    const ds  =  new DataStore(size,{
        columns:columns
    });
    ds.setColumnData("x1",c1);
    ds.setColumnData("x2",c2);
    ds.setColumnData("x3",c3);
    ds.setColumnData("x4",c4);
    return ds;
}

export {getRandomDataStore}