## dataloader function

The asynchronous function takes an array of column objects, the name of the dataSource and the size of the data. it should return an array of objects containing field and data either a SharedArrayBuffer with the  [raw column data]{@tutorial datasource} or a JavaScript array. The latter will be slower as it will have to be converted into native format.
```js
const columnData={
    "x1": [1, 2, 3, 4],
    "color": ["blue", "green", "red", "blue"],
    "tags": ["tag1", "tag1,tag2", "tag3", "tag1,tag3"]
}

async function dataloader(columns, data, size){
    return columns.map(x => ({field: x.field, data: columnData[x.field]}));
}
```

For more complex examples see the [DataLoaders]{@link module:DataLoaders}


## rowDataLoader
This loads arbitrary json data for an item in the datastore and is invoked whenever an item is 'highlighted'. The datastore config should have the `"row_data_loader":true` in order for the loader to be called. An example of a simple rowDataLoader:-
```javascript
async (datasource,index)=>{
    const resp= await fetch(`/my_api?datasource=${datasource}&index=${index}`);
    if (resp.status===200){
        return resp.json();
    }
    return null;
}
```

The data can be accessed via the dataStore's listener 

```javascript
datastore.addListener("my_listener",(type,data)=>{
    if (type==="data_highlighted"){
        //index of the highlighted item
        console.log(data.indexes[0]);
        //data associated with the row
        console.log(data.data)
    }
});
```
Alternatively in a chart override  the `onDataHighlighted()` method:-
```javascript
onDataHighlighted(data){
    this.doSomethingWithData(data.indexes[0],data.data)       
}
```

## binaryDataLoader
Loads binary data given the name of the datasource and the name of the data to load

```javascript
async (datasource,name)=>{
    const resp= await fetch(`/my_api?datasource=${datasource}&name=${name}`);
    if (resp.status===200){
        return resp.arrayBuffer();
    }
    return null;
}
```

The data can then be accessed by the datastore e.g. in a chart
```javascript
async function loadMyData(){
    const data = await this.dataStore.loadBinaryData("my_data");
    //do something with the data
}
```



