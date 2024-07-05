import DataStore from "../datastore/DataStore.js";

const mockColumns= [
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
        field:"colors",
        datatype:"text",
        values:[
        "blue", 
        "purple",
        "red",
        "yellow"
        ],
        colors:[
        "#7FFFD4",
        "#8A2BE2",
        "#DC143C",
        "#FFD700"
        ]
    }
    
];



function getRandomDataStore(size,config={}){
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
            a1[n]=Number.NaN;
        }
        //a1[200]= 420;
    
    }
    const ds  =  new DataStore(size,{
        columns:mockColumns
    });
    ds.setColumnData("x1",c1);
    ds.setColumnData("x2",c2);
    ds.setColumnData("x3",c3);
    ds.setColumnData("colors",c4);
    return ds;
}




function getRandomData(column,size){
    let arrayType=Float32Array;
    let bytes=4;
    if (column.datatype=="unique" || column.datatype=="text"){
        arrayType=Uint8Array;
        bytes=1;
        if(column.datatype=="unique"){
          bytes=column.stringLength;
          arr_len=size*column.stringLength;
        }
       
    }
    const len  = size*bytes;     
    
    const sab = new SharedArrayBuffer(len);
    const new_arr =  new arrayType(sab)
    if (column.datatype === "text"){
      const s = 4;
      for (let n=0;n<len;n++){
          new_arr[n]=Math.floor(Math.random()*s);
      }
    }
    else{
      for (let n=0;n<len;n++){
          new_arr[n]=Math.random()*100;
      }
    }
    return sab

}

function mockDataLoader(){
    return {
        function: (columns, dataSource, size) => {
            return new Promise((resolve,reject)=>{
                const dataList=[];
                for (const column of columns){        
                  dataList.push({data:getRandomData(column,size),field:column.field});
                }       
                resolve(dataList);
            });
        }
    }
}

export function mockDataSource(name, size) {
    return {
        size, name, columns: mockColumns
    }
}

// function mockDataSources(dataStores) {
//     return dataStores.map((ds, i) => {
//         const {size, columns} = ds;
//         return {
//             size,
//             name: "ds" + i,
//             columns
//         }
//     });
// }

export {getRandomDataStore,mockDataLoader,getRandomData, mockColumns}