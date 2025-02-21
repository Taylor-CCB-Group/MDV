# Charts, Views, Links etc

As we develop richer functionality around manipulating Views, we will also want more structured and robust ways of representing them in the code (both Python and JS).

Links are also subject to some further development around API design, as well as UI/UX. Some links are encoded as part of datasources - for example, rows_as_columns defines a particular type of relationship between two datasources, such that matrix data relating to each of them can be manifested. Other links are properties of a particular view - for example, two viv image charts in a view can have their panning/zoom state linked so that as the user interacts with one, the other remains synchronised.


## Chart config objects, lifecycle

We use JSON to describe properties of a chart in a way that can be saved/loaded etc.

At runtime, this JSON object is turned into a mobx `observable`, which means that any React components wrapped in an `observer` should automatically respond to relevant changes - for example, tweaking a slider in a settings dialog will happen in an `action` that sets a new numerical value by mutating that value. The same model can also apply to less primitive properties; assigning new objects to certain properties, altering the contents of arrays, etc. This - in principle - is broadly intended to be the main model of state for charts in future, particularly for any state that should be serialisable (ie, when the user saves, this will be part of the state that will be included - this may also relate to undo in future). Other reactions can be set up in non-react contexts where appropriate, for example in the `func` property of items returned by `chart.getSettings()`.

In practice, however a given chart may model its state internally, what should be universally true, is that the `getConfig()` method should return a JSON object of a form that can be saved to a file, and passed to a chart constructor by `ChartManager` - potentially subject to some validation.

Some charts override `getConfig()` to perform additional serialisation logic, with corresponding deserialisation that applies either in the constructor or elsewhere, while others may operate more directly on `config` (this is considered desirable where it makes sense).

### AddChart dialog

In this initial phase, an internal object with a reduced set of properties is used to represent the aspects of a chart that can be configured here. When it completes, it passes an object in the form of serialised JSON to be processed elsewhere in the system.

The user invokes this in the context of a particular `DataSource` - which dictates what data the chart will have access to, and whether that source is capable of using certain features such as genome browser, image regions etc.

The system looks up entries in `BaseChart.types` (the type of this collection being described in `ChartTypes.ts`) and evaluates which chart types meet that criteria on the basis of any `"required"` specifications, as well as `"allow_user_add"` which can be used to hide things in the front-end. This populates the dropdown used for users to chose the `chartType`.

When a chart type is selected, another set of GUI widgets is populated based on the entries in `chartType.params`, denoting what type (or types) of data is compatible with that entry, and the name that should be shown to the user:

```ts
type Param = "text" | "number" | "multitext" | "text16" | "_multi_column:number" | "_multi_column:all";
...
params?: { type: Param | Param[]; name: string }[]
// e.g. for Dot Plot
[{ type: "text", "name": "Categories" }, { type: "_multi_column:number" }]
```

Note that the names we've denoted `type Param` used for controlling which columns will be shown here do not correspond directly to the strings used for column types. Also, the potential existence of `_multi_column` entries (representing an arbitrary number of values) interspersed into what will ultimately be a flat array, presents a challenge when attempting to process this in an abstract way. In practice, these `_multi_column` entries only appear to exist as the first or last element - still, there is a certain amount of caution is required when altering related code. This is a reason that we may in future favour props having distinct fields (with a revised `"configEntriesUsingFields"` analog) in favour of the current implementation.

As well as these `params`, which appear under a heading "Columns" in the GUI, the system evaluates an optional `chartType.extra_controls` to determine other widgets that should be shown (this is one of a number of places in the codebase where some simple object descriptions are used to describe sets of UI widgets, in similar but not entirely consistent forms - readers without a sufficient disposition towards abstract type generics are liable to be confused by attempting to understand the attempts to formalise these).

The internal state used in this dialog forms a `config` that may be subject to an additional `chartType.init` method, which will mutate to that object, incorporating values from `extra_controls` and transforming it into a more complete form that may be required for the `ChartManager` to be able to deal with it.

In the current implementation of the new `AddChartDialogReact`, following this `init` we perform an additional `serialiseConfig()` - the internal representation prior to this point may have live objects, and the current contract is that the representation understood by `ChartManager` will be the serialised form (this is an implementation detail that could change, but currently while not entirely efficient is at least somewhat coherent).
###### Example of a complex `init()` implementation - Viv charts

From a user-perspective, the only thing to choose when creating a new chart of this type is which `region` should be used, where `region` is a `string` key of `dataSource.regions.all_regions`. This depends on the project being appropriately configured etc such that the chart will be available as an option.

The "internal" representation used in the state of the add chart dialog looks something like this:

```json
{
    "title": "",
    "legend": "",
    "param": [],
    "type": "Viv Scatter Plot (react)",
    "extra": {
        "region": "Unchallenged-SAMPLE_7_ROI_2"
    }
}
```

After calling `init()`, we have

```json
{
    "title": "Unchallenged-SAMPLE_7_ROI_2",
    "legend": "",
    "param": [
        "x",
        "y",
        "leiden"
    ],
    "type": "VivMdvRegionReact",
    "region": "Unchallenged-SAMPLE_7_ROI_2",
    "color_by": "leiden",
    "background_filter": {
        "column": "sample_id",
        "category": "Unchallenged-SAMPLE_7_ROI_2"
    },
    "color_legend": {
        "display": false
    },
    "roi": {
        "min_x": 0,
        "min_y": 0,
        "max_y": 3816,
        "max_x": 3809
    },
    "json": "json/Unchallenged-SAMPLE_7_ROI_2_whole-cell_DAPI_AVGMARKER_200_75.tif.s1.json",
    "background_image": "undefined",
    "radius": 10,
    "viv": {
        "channels": [
            {
                "name": "DAPI"
            }
        ]
    }
}
```

These properties added by `init()` could generally be considered redundant, in the sense that they can be predictably derived from the `"region"`. Many of them are unalterable, others may be changed by the user at runtime (at which point they should certainly be present), and more may also be introduced. There are at present inconsistencies within the way that some of these properties are used internally between different version of the Viv chart which should be reviewed.

It may be preferable to defer the synthesis of such properties to a later stage where possible - in this way we might allow other mechanisms such as the Python API to operate on the simpler initial representation.
### View loading / `ChartManager.addChart()`

Adding a chart (because the "add chart" dialog finished, or while a view is loaded) will call `ChartManager.addChart()`, which internally does some house-keeping such as ensuring that any necessary columns (found in `config.param` as well as `"config[k] for k in chartType.configEntriesUsingColumns"`) are loaded, before calling the constructor for the appropriate chart with a serialised JSON configuration.

Within the chart constructor, a function `initialiseChartConfig(originalConfig, chart)` is called which will deserialise any relevant properties, and perform instrumentation that will allow the chart to respond appropriately to changes.

### User interaction - `chart.getSettings()` etc.

Various arbitrary things may happen during the lifetime of a chart to alter its state - for example, zooming a view might set `config.viewState.zoom`.

One thing that is common to all chart types is some implementation of a method `getSettings()` which returns an array of items to be rendered by a settings dialog. Those can have listeners that execute somewhat arbitrary code, but in many cases they will simply mutate some value in `chart.config`.

For the zoom example, it would be possible for a chart to have both a slider for setting the zoom in the settings dialog, as well as setting it by using relevant gestures on the actual canvas - ideally, the state of all of these should be properly synced (so zooming on the canvas would move the slider).

### Columns, `FieldSpec(s)`, `ColumnQuery`

The architecture previously operated primarily on `string`s which either denoted a column in a datasource by it's `field` property, or a specially structured string like `"Gene scores|GeneName (Gene scores)|0"` which can be interpreted such that it will look up some relevant data and synthesise a column from it.

In order to specify that a given property should respond in an active way to other changes in the system, we use objects in place of these `FieldName` strings. This mechanism may also be used in future to represent other features such as describing columns that will perform some kind of computation - for example synthesising mock-data, or requesting analysis from a server.

So, wherever a chart `config` has a property representing a column, it should be a `FieldSpec`, which may be a `string` or a `ColumnQuery` object. This is a somewhat subtle change, but it is intended to make the system more robust and flexible.

These `ColumnQuery`s must be implemented in a way such that their representation throughout the lifecycle of charts and `config` objects is consistent: there needs to be some way to serialise them (this should be via a `toString()` implementation, such that `JSON.stringify(config)` is sufficient in many cases), and to deserialise them.

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

That means that a chart should be able to accept the deserialised form of a `config` object, and use it to directly access the column data that it represents, in a way that should correctly and transparently handle the loading of column data and reacting to any relevant changes in state (in the implementation, these are `@computed` properties based on `observable` mobx state).

##### `_multi_columns`

It's hard to reason about which entries in `config.param` relate to a given entry - but as far as I've been able to discern, there don't exist any chart types with `_multi_columns` that have more than two entries in `param` - and only one of those is ever `_multi_columns`.

In future, it may be simpler, rather than a `param` array, to have more `"configEntriesUsingColumns"`... so rather than

```js
BaseChart.types["dot_plot"] = {
    "name": "Dot Plot",
    "params": [
        {
            "type": "text",
            "name": "Categories on y-axis"
        },
        {
            "type": "_multi_column:number",
            "name": "Fields on x axis"
        }
    ]
}
```

which means that we need some potentially error-prone logic for reconstituting the columns - in this instance, we know that the first entry is "Categories" and all other entries are "Fields"; but for general-purpose code to look at descriptions of `"params"` in this form and reliably make sense of them is a somewhat awkward task. From observation, we don't have arbitrarily many of these, and any that we do have are always either at the beginning or end of the list... I don't want to contemplate that not being true. We further complicate this with the notion that potentially, multi-column query objects could be interspersed among others within

If we instead had the following (note that this is pseudo-code):

```js
BaseChart.types["dot_plot"] = {
    "name": "Dot Plot",
    "categories": {
		"type": "text",
		"name": "Categories on y-axis"
	},
    "fields": {
		"type": "_multi_column:number",
		"name": "Fields on x axis"
	}
}
```

There would be considerably less ambiguity. The chart itself would be able to refer to `config.categories` safe in the knowledge that it referred to a single column of type `"text"` (`DataColumn<"text">`), and `config.fields` a `DataColumn<NumberDataType>[]`. The more general-purpose code responsible for marshalling a serialised config object into something ready to be digested by a chart (and back into something to be saved) would be spared having to slice-dice (and later glue back together) the flat `config.param` array, which can become wrought with a combinatorial explosion of complexity as the permutations of "is this a string, or an array of strings, or an array with items that might be strings or `DataColumn<T>` where `T` is some generic type derived from the last part of the `"type"` property of the descriptor..." leak into wider.

Old charts are written to expect `"string"` and use it to look up columns in a `DataStore`. New ones prefer to be passed a `DataColumn<>` and forget that `DataStore` even exists - but neither should have to concern themselves too much with where these properties came from.

Note that the above is not precisely the way it should ultimately be written - nor is it a form that would necessarily exactly work in the current code, but it is close. When revisiting how these `BaseChart.types` are specified, we may use `zod` for specifying the shape of these configs in a way that allows for loaded data to be validated and the types to be known by the language-server while editing code (as well as used to generate a JSON-schema which may inform the shaping of the Python API)...


`ColumnQuery` object which can be used to specify a column in a way that it can be resolved to a `Column` object, which will be updated when the underlying data changes. 

The type `FieldSpec` represents what may be a `string` or a `ColumnQuery` - and is used in various places where a column is expected.

#### React

It seems desirable to be able to always assume that where a column is specified in `config`, it will be represented by an actual `DataColumn`; the writer of chart code should not have to litter their logic with `if (typeof param === "string") { ... } else { ... }`. The `"string"` side of that branch entails interaction with `DataStore`, etc...

We have hooks like `useParamColumns()`, which will react to changes in the `config.param`, and return a `DataColumn[]`. The interaction with `DataStore` to load data, reasoning about types are abstracted inside that.

That still leaves other `"configEntriesUsingColumns"` properties that need to be addressed. It may be preferable to have

This may imply that in addition to the deserialisation of `config` currently applied in `BaseChart.constructor`, there could be a further manifestation of deserialised `config` into a more concrete `state`, in which anything originally deriving from a `FieldSpec` - active objects or simpler `FieldName` strings - will be simply present as a `DataColumn`. The chart shouldn't need to worry about the origin of such a thing unless it has reason to.

#### Vanilla JS Charts

These are written in a way that expects columns to be specified as strings representing `column.field` (which might be some specially structured `"Gene scores|GeneName (Gene scores)|0"` which can be interpreted subject to documentation elsewhere).

There was a `"methodsUsingColumns"` in the original architecture which provides a list of methods for which `ChartManager` would instrument extra logic to ensure that when calling these methods, the system would interpret what was passed in the form of `string` would happen only at such a time as the corresponding column-data had been loaded, ready for the chart to use. `"methodsUsingColumns"` can still be used - now superseded by the equivalent with a decorator `@loadColumnData` directly on the method itself.

Both of these methods share a similar implementation, which has now been augmented to react to changes in the `observable` state of active objects passed - so a single call by the user to `DotPlot.setFields([someColumnQueryObject])` will be replaced by something that will call the original `setFields` method with new `FieldName`s in response to changes in relevant state.

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