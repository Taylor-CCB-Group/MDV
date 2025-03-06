import { CTypes, isMultiColumn } from "@/lib/columnTypeHelpers";
import { observer } from "mobx-react-lite";
import { RowsAsColsProps } from "./LinksComponent";
import { DropdownMappedValues, GuiSpec } from "@/charts/charts";
import { useMemo } from "react";
import { action, makeAutoObservable } from "mobx";
import { DropdownAutocompleteComponent } from "./SettingsDialogComponent";
import { g } from "@/lib/utils";
import { getFieldName } from "@/links/link_utils";
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
    const targetColumn = linkedDs.dataStore.columnIndex[name_column];
    if (!targetColumn) {
        throw new Error(`Could not find name_column '${name_column}' in linked dataset`);
    }
    if (linkedDs.size > 2**16) {
        // If we have more than 2^16 rows, I think that means there must be some duplicates...
        // so this is potentially somewhat undefined behaviour.
        console.warn(`Linked dataset has more than 2^16 rows, this may lead to unexpected behaviour`);
    }
    // we need to associate each value with the row index that references it...
    // as far as the code is concerned, this isn't necessarily a 1:1 mapping - but as my understanding of the logic
    // it should be...
    // We may have
    // - more than one row with the same "gene_id"; I think this leads to undefined behaviour, we should give a warning
    // - "gene_id" that doesn't have any corresponding row... maybe this doesn't matter, but we should filter those values out???
    //   (and probably also log a warning)

    // Loop through the values and create the format in pipes '|', by making use of the index.
    //todo - we need to let the user select subgroup...
    //^^ do we allow a mix-match of links/subgroups? We could... but that doesn't mean we should.
    //(we should only do so if we have a clean way of presenting it to the user)
    const sg = Object.keys(subgroups)[0];
    
    // return values in a format that can be used by the dropdown component
    const values = useMemo(() => {
        const labelValueObjects = link.valueToRowIndex.entries().map(([value, rowIndex]) => ({
            label: value,
            fieldName: getFieldName(sg, value, rowIndex)
        }));
        // tricky type here... but 'satisfies' gives us a bit of confidence that we're getting something that will work
        // ! better than `as` typecasting when possible.
        return [[...labelValueObjects], "label", "fieldName"] as const satisfies DropdownMappedValues<"label", "fieldName">;
    }, [targetColumn.values, sg]);
    
    return { values };
}

/**
 * This presents the user with a method for selecting column(s) corresponding to a linked dataset.
 * 
 * There should be a re-usable component for selecting values from a category column, and this should be used here,
 * referring to the `link.name_column` column of the linked dataset.
 * 
 * Values that this component operates on - at the time of writing - will be in a `FieldName` (`string`) format, with pipes separating
 * `subgroup|value (subgroup)|index`, where `subgroup` is the name of the user-selected subgroup, `value (subgroup)` is a human-readable
 * representation of the value, and `index` is the index of the value in the column.
 */
const LinkToColumnComponent = observer(<T extends CTypes,>(props: RowsAsColsProps<T>) => {
    const { linkedDs, link } = props;
    const { name_column, name, subgroups } = link;
    const { values } = useLinkTargetValues(props);
    const { setSelectedColumn, current_value } = props;
    //@ts-expect-error need to review isMultiType logic
    const isMultiType = isMultiColumn(props.type);
    // nb, the <"multidropdown" | "dropdown"> type parameter ends up causing problem vs <"multidropdown"> | <"dropdown">
    // - would be nice to be able to avoid that, this is really difficult to understand and work with.
    // maybe we should have a branch that returns a totally different g() for multidropdown, for example.
    // todo - consider a g<"category"> | g<"multicategory"> for cases such as this where we're chosing row values
    // - this is also used for other things (e.g. density)
    //   (where we don't necessarily need all the faff with associated row indices)
    // - also as well as multi/single chosen options, there is 'multitext' which is a different thing, with a different interface...
    const spec = useMemo(() => makeAutoObservable(g<"multidropdown" | "dropdown">({
        type: 'multidropdown', //!wrong - need to adapt this along with current_value
        // type: isMultiType ? 'multidropdown' : 'dropdown',
        label: `specific '${name}' column`, //todo different label for multiple
        // no ts error! praise be!
        values,
        //this is not what we want to show in a dropdown... this is what a component will be fed if it has opted for 'active selection' mode
        current_value: current_value ? [`${current_value}`] : [],
        func: action((v) => {
            // we don't set current_value here, we use the setSelectedColumn function
            // props.current_value = isMultiType ? v : v[0];
            // const value = `${Object.keys(subgroups)[0]}|${v[0]} (${Object.keys(subgroups)[0]})|${0}`;
            // @ts-expect-error maybe we'll make it so that we have a hook returning a setSelectedLinkedColumn function of narrower type
            setSelectedColumn(isMultiType ? v : v[0])
        })
    })), [values, name, isMultiType, setSelectedColumn, current_value]);
    // const { current_value } = spec;
    return (
        <div className="flex flex-col" style={{ textAlign: 'left' }}>
            {/* really don't want to have this typecast here, what a nuisance! (also not really correct) */}
            <DropdownAutocompleteComponent props={spec as GuiSpec<"dropdown"> | GuiSpec<"multidropdown">} />
        </div>
    )
})

export default LinkToColumnComponent;