## Specifying a View

Views are just collection of charts/widgets for each DataStore and optionally links between the charts. They are are specified in the ChartManager's config

If there is only a single view add the view's config to the ChartManager's config with the key *onlyView*


If there are multiple views, a list of all possible views needs to be defined in the config and a viewloader function in the dataloader config. A viewLoader function is then required in order retrieve a view. This is simply an async function that takes the name of the view and returns the view's config.
```js
async function myViewLoader(viewName){
    return await getView(viewName)
}

const dataLoader={
    function:myDataLoader,
    viewLoader:myViewLoader
}

const config={
    all_views:["view1", "view2", "view3"],
    current_view:"view2",
    columns:myColumns
}

const m = new ChartManager("my-div", myDataSources, dataLoader, config)

```



## The view config

* **dataSources** and object containing the name of each dataSource 
    * layout - either absolute or gridstack (default absolute)
    * panelWidth - the width (as a percentage) of the panel. This is optional,
    as by default each panel will take up an equal amount of space

* **initialCharts** an object containing a list of all the charts for a dataSource

```json
{
    "dataSources":{
        "cells":{
            "layout":"absolute",
            "panelWidth":75
        },
        "genes":{
            "layout":"gridstack",
            "panelWidth":25
        }
    },
    "initialCharts":{
        "cells":[
            {.....},
            {.....}
        ],
        "genes":[
            {......},
            {......}
        ]   
    }
}


```

