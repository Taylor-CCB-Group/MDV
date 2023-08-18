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
