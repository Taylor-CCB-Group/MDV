import { type CTypes, isMultiColumn } from "@/lib/columnTypeHelpers";
import { observer } from "mobx-react-lite";
import type { RowsAsColsProps } from "./LinksComponent";
import type { DropdownMappedValues, GuiSpec } from "@/charts/charts";
import { useMemo } from "react";
import { action, makeAutoObservable } from "mobx";
import { DropdownAutocompleteComponent } from "./SettingsDialogComponent";
import { g, isArray } from "@/lib/utils";
import { getFieldName } from "@/links/link_utils";
/**
 * A hook for managing the state of manually chosen columns from a rows-as-columns link.
 * 
 * Should encapsulate
 * - making sure that an appropriate type of value is selected
 * - providing a way to select the value that has a simple interface, managing the complexity of incoming generic types
 */
function useLinkTargetValues<T extends CTypes, M extends boolean>(props: RowsAsColsProps<T, M>) {
    const { linkedDs, link } = props;
    const { name_column, name, subgroups } = link;
    const targetColumn = linkedDs.dataStore.columnIndex[name_column];
    if (!targetColumn) {
        throw new Error(`Could not find name_column '${name_column}' in linked dataset`);
    }
    if (linkedDs.size > 2**16) {
        // If we have more than 2^16 rows, I think that means there must be some duplicates...
        // so this is potentially somewhat undefined behaviour.
        console.warn("Linked dataset has more than 2^16 rows, this may lead to unexpected behaviour");
    }
    // we need to associate each value with the row index that references it...
    // as far as the code is concerned, this isn't necessarily a 1:1 mapping - but as my understanding of the logic
    // it should be...
    // We may have
    // - more than one row with the same "gene_id"; I think this leads to undefined behaviour, we should give a warning
    // - "gene_id" that doesn't have any corresponding row... maybe this doesn't matter, but we should filter those values out???
    //   (and probably also log a warning)

    // Loop through the values and create the format in pipes '|', by making use of the index.
    //todo - we need to let the user select link/subgroup...
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
        return [[...labelValueObjects], "label", "fieldName"] satisfies DropdownMappedValues<"label", "fieldName">;
    }, [sg, link.valueToRowIndex]);

    const value = useMemo(() => {
        const defaultVal = props.multiple ? [] : values[0][0].fieldName;
        if (!props.current_value) return defaultVal;
        if (props.multiple) {
            // does the current value look like ours? we're not assuming it should:
            // we might not be the active mode...
            if (!isArray(props.current_value)) {
                return defaultVal;
            }
            if (props.current_value.length === 0) {
                return defaultVal;
            }
            const [firstValue] = props.current_value;
            if (typeof firstValue !== 'string') {
                return defaultVal;
            }
            // check that the current value is in the list of possible values
            if (!firstValue.includes("|") || !link.valueToRowIndex.has(firstValue)) {
                return defaultVal;
            }
            return props.current_value;
        } else {
            // if it's an array, we're not the active mode
            if (isArray(props.current_value)) {
                return defaultVal;
            }
            if (typeof props.current_value !== 'string') {
                return defaultVal;
            }
            // check that the current value is in the list of possible values
            if (!props.current_value.includes("|") || !link.valueToRowIndex.has(props.current_value)) {
                return defaultVal;
            }
            return props.current_value;
        }
    }, [props.current_value, link.valueToRowIndex, props.multiple, values[0]]);

    const { setSelectedColumn } = props;
    const isMultiType = props.multiple;
    
    
    
    return { values, value };
}

/**
 * This presents the user with a method for selecting column(s) corresponding to a linked dataset.
 * 
 * There should be a re-usable component for selecting values from a category column, and this should be used here,
 * referring to the `link.name_column` column of the linked dataset... at the moment, we deal with this in a more ad-hoc way.
 * 
 * Values that this component operates on - at the time of writing - will be in a `FieldName` (`string`) format, with pipes separating
 * `subgroup|value (subgroup)|index`, where `subgroup` is the name of the user-selected subgroup, `value (subgroup)` is a human-readable
 * representation of the value, and `index` is the index of the value in the column.
 */
const LinkToColumnComponent = observer(<T extends CTypes, M extends boolean>(props: RowsAsColsProps<T, M>) => {
    const { linkedDs, link } = props;
    const { name_column, name, subgroups } = link;
    const { values, value } = useLinkTargetValues(props);
    const { setSelectedColumn } = props;
    const isMultiType = props.multiple;
    
    // nb, the <"multidropdown" | "dropdown"> type parameter ends up causing problem vs <"multidropdown"> | <"dropdown">
    // - would be nice to be able to avoid that, this is really difficult to understand and work with.
    // maybe we should have a branch that returns a totally different g() for multidropdown, for example.
    // todo - consider a g<"category"> | g<"multicategory"> for cases such as this where we're chosing row values
    // - this is also used for other things (e.g. density)
    //   (where we don't necessarily need all the faff with associated row indices)
    // - also as well as multi/single chosen options, there is 'multitext' which is a different thing, with a different interface...
    
    /// we don't want a new spec every time the value changes... we want to give it an observable value
    // (I think)- definition of this spec may want to be in the hook, with an observable value that isn't
    // the real props.current_value???
    const spec = useMemo(() => makeAutoObservable(g<"multidropdown" | "dropdown">({
        type: isMultiType ? 'multidropdown' : 'dropdown',
        label: `specific '${name}' column`, //todo different label for multiple
        // no ts error! praise be!
        values,
        current_value: value,
        func: action((v) => {
            // @ts-expect-error maybe we'll make it so that we have a hook returning a setSelectedLinkedColumn function of narrower type
            setSelectedColumn(v);
        })
    })), [values, name, isMultiType, setSelectedColumn, value]);
    // const { current_value } = spec;
    return (
        <div className="flex flex-col" style={{ textAlign: 'left' }}>
            {/* really don't want to have this typecast here, what a nuisance! (also not really correct) */}
            <DropdownAutocompleteComponent props={spec as GuiSpec<"dropdown"> | GuiSpec<"multidropdown">} />
        </div>
    )
})

export default LinkToColumnComponent;