# Charts, Views, Links etc

As we develop richer functionality around manipulating Views, we will also want more structured and robust ways of representing them in the code (both Python and JS).

Links are also subject to some further development around API design, as well as UI/UX. Some links are encoded as part of datasources - for example, rows_as_columns defines a particular type of relationship between two datasources, such that matrix data relating to each of them can be manifested. Other links are properties of a particular view - for example, two viv image charts in a view can have their panning/zoom state linked so that as the user interacts with one, the other remains synchronised.


## JavaScript / TypeScript

### Chart config objects, lifecycle

We use JSON to describe properties of a chart in a way that can be saved/loaded etc.

At runtime, this JSON object is turned into a mobx `observable`, which means that any React components wrapped in an `observer` should automatically respond to relevant changes in config - for example, tweaking a slider in a settings dialog will happen in an `action` that sets a new numerical value. This - in principle - is broadly intended to be the main model of state for charts, particularly for any state that should be serializable (ie, when the user saves, this will be part of the state that will be included - this may also relate to undo in future). Other reactions can be set up in non-react contexts where appropriate, for example in the `func` property of items returned by `chart.getSettings()`.

In practice, however a given chart may model its state internally, what should be universally true, is that the `getConfig()` method should return a JSON object of this form, and that is the form passed to the chart constructor by `ChartManager` - potentially subject to some validation.

Some charts override `getConfig()` to perform additional serialisation logic, with corresponding deserialisation that applies either in the constructor or elsewhere, while others may operate more directly on `config` (this is considered desirable where it makes sense).


### Columns, `FieldSpec(s)`, `ColumnQuery`

The architecture previously operated primarily on `string`s which either denoted a column in a datasource by it's `field` property, or a specially structured string like `"Gene scores|GeneName (Gene scores)|0"` which can be interpreted such that it will look up some relevant data and synthesise a column from it.

In order to specify that a given property should respond in an active way to other changes in the system, we use objects in place of these `FieldName` strings. This mechanism may also be used in future to represent other features such as describing columns that will perform some kind of computation - for example synthesising mock-data, or requesting analysis from a server.

So, wherever a chart `config` has a property representing a column, it should be a `FieldSpec`, which may be a `string` or a `ColumnQuery` object. This is a somewhat subtle change, but it is intended to make the system more robust and flexible.

These must be implemented in a way such that their representation throughout the lifecycle of charts and `config` objects is consistent: there needs to be some way to serialise them (this should be via a `toString()` implementation, such that `JSON.stringify(config)` is sufficient in many cases), and to deserialise them.

So, the serialised form may look something like this:

```json
{
  "type": "RowsAsColsQuery",
  "linkedDsName": "Genes",
  "subgroup": "Gene scores",
  "maxItems": 10
}
```

This will be deserialised into a `MultiColumnQuery` object, with the following properties:

```typescript
interface MultiColumnQuery {
    columns: DataColumn<DataType>[]; //could have a generic type for this?
    fields: FieldName[]; //`FieldName` is just an alias for `string`, but denotes that it refers to a column
}
```

That means that a chart can accept the deserialised form of a `config` object, and use it to directly access the column data that it represents, in a way that should correctly and transparently handle the loading of column data and reacting to any relevant changes in state (in the implementation, these are `@computed` properties based on `observable` mobx state).

Charts which are not written in a way that they can recognises these objects can be instrumented to do so - meaning that they still internally accept `string`.


`ColumnQuery` object which can be used to specify a column in a way that it can be resolved to a `Column` object, which will be updated when the underlying data changes. 

The type `FieldSpec` represents what may be a `string` or a `ColumnQuery` - and is used in various places where a column is expected.

#### React



#### Vanilla JS Charts

These are written in a way that expects columns to be specified as strings representing `column.field` (which might be some specially structured `"Gene scores|GeneName (Gene scores)|0"` which can be interpreted subject to documentation elsewhere).

There was a `"methodsUsingColumns"` in the original architecture which provides a list of methods for which `ChartManager` would instrument extra 

#### VivMDVReact

This uses a set of hooks adapted from the original Avivator code-base, with Zustand for state management wrapped in an additional layer of React’s context API such that each chart can have it’s own model of that state.

This is then somewhat marshalled into config.viv in getConfig() - and there are also ways in which linked charts can access and manipulate some aspects of this state (subject to review). Note that the structure of this is different to the older VivScatterPlot - and not necessarily ‘better’ (this should be the subject of a separate, more in-depth document).

GUI Widgets - add chart vs getSettings() etc.

## Regions

When the user creates something like a spatial chart, they don’t manually choose param values - those are handled 

## Wider state management

In future, we may consider wider review of the architecture for state management, such that we can support features like undo/redo in a consistent way, or allow for dynamically updating the state of datasources metadata, for example (as may happen if new data is added on the server at runtime, and we want to view it without refreshing the page).

We may move toward wider rendering of components with React - the notion that ‘interface is a function of state’ may mean that the role of ChartManager, and indeed classes derived from BaseChart, become less central to the architecture.

If currentView is an observable store, it ought to be possible that a functional-react component responsible for rendering the entire application DOM would appropriately respond to changes in its structure - such as a chart being added, deleted, moved, etc.

# Python

The current implementation of Chart APIs follows a pattern, with a BasePlot inherited by other types of charts. It is mostly a wrapper for the plot_data: dict. 


```python
class BasePlot:
    def __init__(self, title, plot_type, params, size, position, id=None, **kwargs):
        self.plot_data = {
            "title": title,
            "type": plot_type,
            "param": params,
            "size": size,
            "position": position,
            "id": id if id else self.generate_id(),
        }
        # arbitary key-value pairs can be added to the plot_data
        for key, value in kwargs.items():
            self.plot_data[key] = value
    def generate_id(self):
        """Generate a unique ID for the plot"""
        return str("".join(random.choices(string.ascii_letters, k=6)))
    def set_legend(self, legend):
        self.plot_data["legend"] = legend

class DotPlot(BasePlot):
    def __init__(self, title, params, size, position, id=None):
        super().__init__(title, "dot_plot", params, size, position, id)
    def set_axis_properties(self, axis, properties):
        if "axis" not in self.plot_data:
            self.plot_data["axis"] = {}
        self.plot_data["axis"][axis] = properties
    def set_color_scale(self, log_scale=False):
        self.plot_data["color_scale"] = {"log": log_scale}
    def set_color_legend(self, display, position):
        self.plot_data["color_legend"] = {"display": display, "pos": position}
    def set_fraction_legend(self, display, position):
        self.plot_data["fraction_legend"] = {"display": display, "pos": position}
    # Additional methods for customization (e.g., tooltip visibility) can be added here
```
As far as I’m aware, other then setting values in parts of that structure, the only logic that it performs at this stage is automatically generating an id when not provided.

To use it, you write something like this:


```python
# Get the cells dataframe
cell_df = p.get_datasource_as_dataframe("cells")
# Add a table for displaying metadata
table_plot = TablePlot(
    title="Metadata Table",
    params=list(cell_df.columns),
    size=[600, 500],
    position=[850, 10],
)
# Configure the project view with the table
view_config = {
    "initialCharts": {
        "cells": [
            table_plot.plot_data
        ]
    }
}
p.set_view("default", view_config)
```

In order to maintain the API, we edit the corresponding classes manually - it may be that we could automate that from a JSON-API or similar spec at some point.

The API could also be made more concise using TypedDict (as long as we don’t want to introduce any extra internal logic such as validation), or autoclass (used in Vitessce) etc.