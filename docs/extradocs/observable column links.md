In the file `link_utils.ts`, there is an exported method `getRowsAsColumnsLinks(dataStore)`.

This will add an `observableFields` property to `link` objects held by a `DataStore` - with listeners that respond to `"data_highlighted"` and `"filtered"` events from the linked `DataStore` and update that array with values from the selected rows in a `mobx action`.

The interface `IRowAsColumn` defining the type of entries in `observableFields` defines that these will have the following properties:

```ts
interface IRowAsColumn {
    index: number;
    name: ColumnName;
    fieldName: FieldName;
    column: DataColumn<DataType>;
}
```

Of these, only the `index` is initially used by the `RAColumn` constructor, with the others being `@computed` properties that are evaluated lazily based on that index when required.

The `observableFields` array can be arbitrarily large - the code that consumes these values in order to know which columns are of interest will generally specify how many items it is able to represent. So while there may be many instances of `RAColumn`, it is only when the `@computed` values `fieldName` and `column` are accessed that they have significant cost. These properties encapsulate the logic for formatting `FieldName` in the format expected by `DataStore.addColumnFromField(fieldName)`, and calling that method to manifest the appropriate `DataColumn` object. The actual data loading is not implicitly handled by this class, but should be exposed in a way that will be easy to do given a reference to a `DataColumn` object without needing to directly interface with `DataLoader` / `DataStore` / `ChartManager` etc.

We may later consider making `observableFields` itself an `iterator` or similar, such that rather than potentially having to process a large number of items many of which won't be used, it could lazily evaluate the needed options on demand. However, the objects created are very tiny