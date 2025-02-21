# Approach

The original code is written in JavaScript, with some `jsdoc` annotations that can provide some level of hint to a language server.

Gradually, there has been a movement to formalise these into TypeScript. At times these types - in their effort to express somewhat dynamic generic variance etc - have become undeniably complex, and may not be easy to immediately grasp from reading the code that defines them. It is also true that the organisation of this code follows a somewhat haphazard structure.

Even where that is true (and it is hoped that we can address this), the aim is that the code which uses them is able to easily **infer** correctly how values passed around are shaped, what they represent and so on, so that the author of an individual component should have rich contextually relevant information, and be able to confidently use the properties of these objects safe in the knowledge that the responsibility for checking their conformance to relevant protocols will be appropriately delegated, and that potentially unsafe or 'incorrect' operations will be flagged.

They should be able to "go to definition" to follow these relations through, and hopefully find some helpful doc-strings along the way.

We aim to trap errors at the earliest possible time, so that if - for instance - a chart is loaded with a configuration that attempts to assign a categorical column where it is required to be numeric, the user should be presented with an error to that effect - rather than one which arises at some later stage, when the chart attempts to use the `minMax` stats from the column metadata only to find that property is `undefined` - raising an exception that bears little coherent relation to the real cause. Types can help us to ensure that large categories of error are caught at the time of code-editing, although this does need to be weighed against .

Of course, TypeScript itself can only decorate the code with assertions that we make as programmers - but it should be able to give us the facility to express what we deem to be knowable about the state of the system at some time: perhaps we are sure that the "add chart" dialog will never produce a configuration that would lead to an error such as that described above, but if we are loading view data from a server, we would like to reason about the fact that we cannot be certain of it's correctness until some form of validation occurs - and this reasoning can be expressed by passing some value of a type which has known and verified characteristics. In this way, as the assertions encoded in these characteristics change through the introduction of new features or refactoring, we can track which parts of the code will be impacted - to the extent that our use of types is correct.

This can avoid the author of a chart being encumbered with checking - for instance - whether a value representing "a column I need to render" is an array, or a FieldName (string) which can be used via reference a DataStore, whether the data needs to be fetched from the server, has a compatible data type - or with more recent additions, how to know when the chosen option has changed, etc.

## Config

Charts generally extend `BaseChart<T extends BaseConfig>` where `BaseConfig` has some common props: `id: string, size: [x: number, y: number], type: string, ...`, and other properties which can be composed for a particular type of chart implementation.

The `BaseChart.constructor(dataStore: DataStore, div: HTMLDivElement | string, config: T)` takes an object `T` of the kind of config that this chart will understand... however, what this definition doesn't currently clearly express is the contract we've stated elsewhere that there is a distinction between the **serialised form** of a given config, vs the **live state** used by an instance of a chart, and further to that, the **concrete resolution** of which columns (or perhaps other derived state) are the result of the specified configuration at a given moment. That latter part (referred to here as **concrete resolution** of state) is further delineated based on whether the associated data is loaded yet, and whether the rest of the code for a given chart is expecting to have these manifest as `FieldName (string)` or `DataColumn`. Perhaps it is useful that these latter manifestations should be referred to as `state` as differentiated from `config` and in turn `serialisedConfig`. Currently, existing chart implementations tend to refer to `this.config` internally to access relevant state, and there is a reasonably large surface-area to consider in order to allow them to seamlessly adapt to a different version of this system.

They should be able to operate as much as possible on the "concrete" level, with the responsibility for aspects such as how `reactions` might be wired to a `RowsAsColsQuery` instance being abstracted away.

So, what we aim to clarify is the way we express this categorisation, and the protocol by which we synthesise the appropriate state, such that the system can operate as seamlessly and reliably as possible, with as little burden placed on a particular chart implementation as possible, and with minimal disruption to existing code.


In general, I would argue that it is more convenient to write

```ts
const xData = config.x.data; // where x can be inferred to be a LoadedDataColumn<NumberDataType>
```

or perhaps

```ts
await config.x.data; // we could have lazily-loaded column data when a sub-chart comes into view etc.
// in this case, the `data: Promise<Float32Array>` should have a getter that is responsible for dealing with the DataStore.
// this should be a minor change to the existing implementation
```

vs

```ts
const x = config.param[0]; // it's easy-ish to know that param[0] is 'x' here... but maybe sometimes we need to refer back to the way the BaseChart.type is specified to check. It's less easy in general with _multi_column etc...
// we can probably find a reference to a dataStore... but I'd like not to have to think about it here.
// at least we can hopefully assume that whatever brought us here will ensure the data is loaded at this point
const xData = dataStore.columnIndex[x].data; //oops, what if x is some "gs|TGF1 (gs)|0" kinda value... that might be ok (I think?) as long as something else in the system has loaded it so it'll be available in columnIndex... if we find ourselves getting to here with x not being a string, though, then it's game-over.
```

This example somewhat conflates a few concerns: 
- convention of passing string (or perhaps some representation of a live-query) vs column object - and how the underlying type for that should be expressed so that a given chart can be written to avoid falling back on a certain class of "`if (typeof x === "string") ...` " statements (which may seem simple in a particular individual case, but make reasoning about the impact of changes to the system hard, and can somewhat explode when the permutations of `string | string[] | ...` grows).
- using `config.param` with loosely held associations to their meanings vs more semantically named properties - which is something for later review.
- where the responsibility for the data-loading is held, in situations where various permutations of active state may be involved.

With React-based code using `param`, we often currently have something like this:

```ts
const xData = useParamColumns()[0].data; //this should abstract away everything apart from needing to know that we want the data from param[0], including reaction to relevant changes.
```

For illustration, pseudo-code will be used that is less based on the `param` array, although in practice, as of this writing most input data for charts tends to be defined through `params` in the corresponding `BaseChart.types` entry - which is a kind of schema defining how the chart will be represented in "Add Chart Dialog", such as

```ts
BaseChart.types["wgl_scatter_plot"] = {
	"name": "2D Scatter Plot",
	"params": [
		{
			"type": "number",
			"name": "X axis"
		},
		{
			"type": "number",
			"name": "Y axis"
		}
	]
}
```


Ideally, a config `T` describing the shape and type of data for a given chart type be could defined once, and used as a generic parameter of another type representing a given lifecycle phase (preferably such that the associated generic types are both human and machine readable). It would be desirable for the seemingly somewhat arbitrary distinction between which things are configurable at the "add chart" stage, vs through `getSettings()`, was also based on this single-source-of-truth... Not to mention the Python API.

The concern about what properties are able to be changed during the lifetime of a chart as opposed to when it is constructed are brought to the fore with the introduction of values that are able to intrinsically vary.

Suppose we have something like this:

```ts
// note - this is pseudocode, related but not identical to how things are currently expressed
// in place of the TypeScript type here - which decomposes to nothing at runtime - it may be written as a zod schema (or in some other way that can manifest as code that will facilitate validation, transformation etc)
type ScatterPlotConfig = {
  x: ColumnSpec<number>;
  y: ColumnSpec<number>;
  radius: number;
  opacity: number;
  ...
} & BaseConfig
```

This may manifest in the system as something like this:

```ts
const configWithFieldNames = {
  x: "x", y: "y", ... //where "x" & "y" correspond to some column.field in the dataSource
} satisfies ConcreteFieldConfig<ScatterPlotConfig>;
```

or perhaps an `observable`

```ts
const config: ConcreteDataColumnConfig<ScatterPlotConfig> = {
	x: LoadedDataColumn<NumberDataType>,
	y: LoadedDataColumn<NumberDataType>,
	...
}
```

Also - a config like this

## Columns

A fundamental aspect of the system is columnar data and how it is represented, with associated metadata.

