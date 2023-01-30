## specifying a view

Views are just collection of charts/widgets for each DataStore and optionally links between the charts. They are are specified in the ChartManager's config

If there is only a single view then initialViews


If there are multiple views, a list of all possible views needs to be defined in the config and a viewloader function in the dataloader config. Optionally and initial 



```
{
    "datasource1":{
        "initialCharts":[

        ],
        "links":[

        ]
    },
    "datasource2":{
        "initialCharts":[

        ],
        "links":[

        ]
    }
}


```

If neither of the above are in the config then a blank view for each DataStore will be shown



###

an object containing a key 

 If you do not specify any views then the app will just show a blank canvas for each dataset
If there is only a single view then it can be specified in the config with onlyView.
Otherwise, you need to specify the available views and the current view in the config

```

    {
        "all_views":[
            "view 1",
            "view 2",
            "view 3"
        ],
        "initial_view":"view 1"
    }
```
A viewLoader is then required to load in the views simply an async function that takes the name of the view and returns the views config
```
async function viewLoader(viewName){
    return await getView(viewName)
}

```