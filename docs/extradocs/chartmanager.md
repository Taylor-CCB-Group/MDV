
The ChartManager object links DataSources, loads column data on demand and manages charts/widgets. It is the main interface for interacting with MDV. The ChartManager's constructor requires list of DataStore configs, a dataloader object, a config for construction and optionally a listener. 
```
    const cm =  new ChartManager("mydiv", datasources, dataloader, config, listener)
```
* **mydiv** - the id or the element to house the app
* **datasources** - list of DataStore [configs]{@tutorial datasource}, these configs should also include a size parameter showing the number of rows in the data set
* **dataloader** - a dataloader which comprises of three parameters
    * **function** - this [function]{@tutorial dataloader} accepts a list of columns and returns a promise (not needed if all data is to be loaded from files)
    * **viewLoader** - a function that will return the [view]{@tutorial views} given the view name
    * **rowDataLoader** - an optional function that accepts a datasource name and index and returns an an object with data. The
    function will be called if the datasource's config contains row_data_loader:true and DataHighlighted is invoked  
    * **files** - specifies static files (tsv, csv or json), which contain the data to display. Useful for small amounts of data (100 000s rows) and testing. If all the data is to be loaded dynamically then this is not required.
    ```json
    [
        {
            "type":"tsv",
            "dataSource":"cells",
            "url":"data/cell_all_archr.tsv"
        },
        {
            "type":"tsv",
            "dataSource":"genes",
            "url2":"data/genes.txt"
        }
    ]
    ```
* **config** - A config 
* **listener** An optional listener function, although this can be added later with the *addListener* method.



## Listeners

Can be added with the method *addListener(id, function)* and removed with *removeListener(id)*. Alternatively a listener can be added as the last parameter when constructing the ChartManager object.

The listener should be a function which receives the type of event, the ChartManager object and any data associated with the event. A typical listener would be:-

```js
    (type, cm, data)=>{
        switch(type){
            case "view_loaded":
                ..do stuff with data
                break;
            case "state_saved":
                ..push data to server
                break;
        }
    }
```

The types of listeners are:- 

* **chart_added**  Called When a chart is added with notify=true e.g. when a user adds a chart. The data received is the chart object itself.
* **chart_removed** Called when a chart is removed with notify=true e.g. when a user removes a chart. The data received is the chart object itself.
* **state_saved** Called when the user saves the state. The data being the state object
* **view_loaded** Called when a view has been completely loaded i.e. all data retrieved and all the charts added. The data being passed is the view that was loaded
* **filtered** Called when a DataStore is filtered, passing the Dimension that has done the filtering



## state_saved

Called when the user saves the data. The object passed to the listener consists of the following:-

* view - A config containing all data for the view
* currentView - 
* all_views -  A list of views 
* updatedColumns -  a dictionary with a entry for each datasource 
    * columns - a list of all columns that have either been updated or added containing:
        * metadata- the columns metadata
        * data an array of the column's raw data
    * added -  a list of all column id(fields) that have been added
    * removed -  a list of all columns that have been deleted
    * colors_changed -a list of columns whose color scheme has changed
        * column - the id of the column
        * colors -the new color scheme (list of hex values)

* metadata - a dictionary of all datasources whose metadata has been updated, with each entry being a dictionary of the parameter and value

```json
{
    "view":{....},
    "cuurentView":"default",
    "all_views":["default","myview"],
    "updatedColumns":{
        "cells":{
            "columns":[
                {
                    "metadata":{"field":"new annotations",...},
                    "data":[2, 1, 2.....]
                },
                {
                    "metadata":{"field":"annotations",...},
                    "data":[0, 0, 2.....]
                }
            ],
            "added":["new annotations"],
            "removed":[],
            "colors_changed":[
                {
                    "column":"annotations",
                    "colors":["#4532FF",....]
                }
            ]
        }
    },
    "metadata":{
        "cells":{
            "param1":{....},
            "param2":{....}
        }
    }

}


```