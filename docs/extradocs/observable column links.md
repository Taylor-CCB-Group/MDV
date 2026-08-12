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

## What these links are, and where they apply

**Any numeric column parameter should be able to accept one of these links.** That is the rule the
column picker implements: `ColumnSelectionComponent` offers the "active link" tab wherever
`paramAcceptsNumeric` holds, and a consumer that offers the tab is expected to handle a
`RowsAsColsQuery` coming back rather than only a column name.

As currently used, a `RowsAsColsLink` is really **"choose a `var` from a table"** — the linked
datasource is a gene/feature table and the query picks a column of the expression matrix. Reading
it that way is worth keeping in mind for two things that do not exist yet: columns computed from an
expression graph (the same problem — a column identity that resolves late and can change while a
chart is open), and fetching a var for an annotated element that covers **far fewer rows than the
datasource**, where materialising the whole column is mostly waste. See
[docs/design/spatial-tables/00-table-element-association.md](../design/spatial-tables/00-table-element-association.md#columns-that-come-from-elsewhere-links-vars-computed).