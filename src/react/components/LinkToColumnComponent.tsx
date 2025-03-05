import { CTypes, isMultiColumn } from "@/lib/columnTypeHelpers";
import { observer } from "mobx-react-lite";
import { RowsAsColsProps } from "./LinksComponent";
import { DataColumn, DataType, GuiSpec } from "@/charts/charts";
import { useMemo } from "react";
import { action, makeAutoObservable } from "mobx";
import { DropdownAutocompleteComponent } from "./SettingsDialogComponent";
import { g } from "@/lib/utils";
/**
 * A hook for managing the state of manually chosen columns from a rows-as-columns link.
 * 
 * Should encapsulate
 * - making sure that an appropriate type of value is selected
 * - providing a way to select the value that has a simple interface, managing the complexity of incoming generic types
 */
function useLinkTargetValues<T extends CTypes,>(props: RowsAsColsProps<T>) {
    const { linkedDs, link } = props;
    const { name_column, name, subgroups } = link;
    const cm = window.mdv.chartManager;
    const targetColumn = cm.getDataSource(linkedDs.name).columnIndex[name_column] as DataColumn<DataType>;
    // we need to associate each value with the row index that references it...
    // as far as the code is concerned, this isn't necessarily a 1:1 mapping - but as my understanding of the logic
    // it should be...
    // We may have
    // - more than one row with the same "gene_id"; I think this leads to undefined behaviour, we should give a warning
    // - "gene_id" that doesn't have any corresponding row... maybe this doesn't matter, but we should filter those values out???
    //   (and probably also log a warning)
    // - If we have more than 2^16 rows, I think that means there must be some duplicates...

    // Loop through the values and create the format in pipes '|', by making use of the index.
    //todo - what if we have multiple subgroups, need to take care of that the target column is not returning subgroups
    const formattedValues = targetColumn.values.map((value, index) => (`${Object.keys(subgroups)[0]}|${value} (${Object.keys(subgroups)[0]})|${index}`))
    
    return {values: targetColumn.values, formattedValues}; //this output should be in the 'FieldName' form with pipes
}

const LinkToColumnComponent = observer(<T extends CTypes,>(props: RowsAsColsProps<T>) => {
    const { linkedDs, link } = props;
    const { name_column, name, subgroups } = link;
    const {values, formattedValues} = useLinkTargetValues(props);
    const { setSelectedColumn, current_value } = props;
    //@ts-expect-error need to review isMultiType logic
    const isMultiType = isMultiColumn(props.type);

    // nb, the <"multidropdown" | "dropdown"> type parameter ends up causing problem vs <"multidropdown"> | <"dropdown">
    // - would be nice to be able to avoid that, this is really difficult to understand and work with.
    // maybe we should have a branch that returns a totally different g() for multidropdown, for example.
    const spec = useMemo(() => makeAutoObservable(g<"multidropdown" | "dropdown">({
        // type: isMultiType ? 'multidropdown' : 'dropdown',
        type: 'multidropdown',
        // name: name_column,
        label: `specific '${name}' column`, //todo different label for multiple
        values: [values],
        //this is not what we want to show in a dropdown... this is what a component will be fed if it has opted for 'active selection' mode
        current_value: [values[0]],
        func: action((v) => {
            // we don't set current_value here, we use the setSelectedColumn function
            // props.current_value = isMultiType ? v : v[0];
            // const value = `${Object.keys(subgroups)[0]}|${v[0]} (${Object.keys(subgroups)[0]})|${0}`;
            // @ts-expect-error maybe we'll make it so that we have a hook returning a setSelectedLinkedColumn function of narrower type
            setSelectedColumn(isMultiType ? v : v[0])
        })
    })), [values, name, isMultiType, setSelectedColumn]);
    // const { current_value } = spec;
    return (
        <div className="flex flex-col" style={{ textAlign: 'left' }}>
            {/* really don't want to have this typecast here, what a nuisance! (also not really correct) */}
            <DropdownAutocompleteComponent props={spec as GuiSpec<"dropdown"> | GuiSpec<"multidropdown">} />
        </div>
    )
})

export default LinkToColumnComponent;