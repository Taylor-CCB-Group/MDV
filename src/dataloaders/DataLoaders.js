/**
*  A collection of dataloaders and helper functions that can
* can be used with {@link ChartManager}
* @module DataLoaders
*
* 
**/


/**
* @memberof module:DataLoaders
* @param {ArrayBuffer} - an array buffer containing raw concatenated
* column data
* @param {object} columns - any array of column objects 
* <ul>
*   <li> field </li>
*   <li> datatype </li>
*   <li> sgtype </li>
* </ul>
* @param {integer} - the size of the columns
* @returns {object[]} a list of objects containing each colums's field name
* and data
**/
function processArrayBuffer(data,columns,size){
    const dataList= [];
	let offset=0;
	for (let column of columns){
        //default values for numbers
		let arrayType= column.datatype==="int32"?Int32Array:Float32Array;
        //the length of the typed array
		let arr_len=size;
        //the number of bytes for each item in the column's data
 		let bytes=4;
        //set the correct values for other datatypes
		if (column.datatype=="unique" || column.datatype=="text"){
			arrayType=Uint8Array;
 			bytes=1;
			if(column.datatype=="unique"){
				bytes=column.stringLength;
				arr_len=size*bytes;
			}
		}
        else if (column.datatype==="multitext"){
            arrayType=Uint16Array;
            bytes=column.stringLength*2;
            arr_len=size*column.stringLength;
        }
        //special way to deal with sparse data 
        //assumes all sparse data is integer/double
        //the first 4 bytes specifies the number of values (n)
        //The next n*4 bytes are the indexes of these values (Uint32)
        //the next n*4 bytes are the actual values (Float32) 
		if (column.sgtype=="sparse"){
			//first byte is length of data
			const l = new Uint32Array(data,offset,1)[0];
			offset+=4;
            //get the indexes
			const indexes = new Uint32Array(data,offset,l);
			offset+=l*4;
            //get the values
			const values = new Float32Array(data,offset,l)
			offset+=l*4;
			const sab = new SharedArrayBuffer(size*4);
			const new_arr= new Float32Array(sab);
            //fill array with missing values
			new_arr.fill(NaN);
            //add the sparse data
			for (let i=0;i<indexes.length;i++){
				new_arr[indexes[i]]=values[i];	
			}
			dataList.push({data:sab,field:column.field});
		}
		else{
			const len  = size*bytes;
            //get the data from the arraybuffer into a SharedArrayBuffer
            //unfortunately cannot be copied directly-  have to via a typed Array
			let arr = new arrayType(data,offset,arr_len);
			const sab = new SharedArrayBuffer(len);
			const new_arr =  new arrayType(sab)
			new_arr.set(arr,0);
			dataList.push({data:sab,field:column.field})
			offset+=len;
		}		
	}
	return dataList;
}




/**
* Gets bytes from an API. The data loader will send a post request 
* to the url with with a jsonified object containing the datasource
* and column information 
* <pre>
* {
*    "data_source":"mydataource"
*    "columns":[{"field":"x1","datatype":"integer"}]
* }
* </pre>
* returns a d
* @memberof module:DataLoaders
* @param {string} - The url of the api
* @returns {function} a dataloader that can be used to construct {@link ChartManager}
**/
function getArrayBufferDataLoader(url){
	return async(columns,dataSource,size)=>{
		//get the data
		const response = await fetch(url,{
			method:"POST",
			body:JSON.stringify({columns:columns,data_source:dataSource}),
			headers: {
				'Content-Type': 'application/json'
			}
		});
		//the data is any arraybuffer containing each individual 
		//column's raw data
		const data = await response.arrayBuffer();
		return processArrayBuffer(data,columns,size)
	}
}


export {getArrayBufferDataLoader,processArrayBuffer};