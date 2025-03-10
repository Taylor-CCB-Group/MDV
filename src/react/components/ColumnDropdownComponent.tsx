import { columnMatchesType, inferGenericColumnGuiProps } from "@/lib/columnTypeHelpers";
import type { ColumnSelectionProps, CTypes } from "@/lib/columnTypeHelpers";
import { useDataStore } from "../context";
import type { DataColumn, DataType } from "@/charts/charts";
import { useMemo } from "react";
import { observer } from "mobx-react-lite";
import { Autocomplete, Checkbox } from "@mui/material";
import { isArray } from "@/lib/utils";
import { TextFieldExtended } from "./TextFieldExtended";
import Grid from '@mui/material/Grid2';
import CheckBoxOutlineBlankIcon from "@mui/icons-material/CheckBoxOutlineBlank";
import CheckBoxIcon from "@mui/icons-material/CheckBox";

const icon = <CheckBoxOutlineBlankIcon fontSize="small" />;
const checkedIcon = <CheckBoxIcon fontSize="small" />;

// todo - all kinds of things to do here...
// - fix multi columns
// - column groups

/**
 * Check if a column selection prop is for multiple columns...
 */
function isMultiColProp<T extends CTypes, M extends boolean>(p: ColumnSelectionProps<T, M>): M {
    return p.multiple;
}

const useColumnDropdownValue = <T extends CTypes, M extends boolean>(gProps: ColumnSelectionProps<T, M>) => {

    const props = inferGenericColumnGuiProps(gProps);
    const { setSelectedColumn, placeholder, type, current_value } = props;
    
    const isMultiType = isMultiColProp(props);
    // are we sure that props.dataStore will be right 
    // - are there any inconsistencies between this and the current DataStoreContext?
    const dataStore = useDataStore(props.dataStore);
    // todo column groups
    const columns: DataColumn<DataType>[] = useMemo(
        () => dataStore.columns
            .filter((c) => !props.exclude?.includes(c.name))
            .filter((c) => !c.field.includes("|")) //exclude linked columns the hacky way for now
            .filter((c) => columnMatchesType(c, type))
            ,
        [dataStore, props.exclude, type],
    );

    // todo allow "None" if optional

    const value = useMemo(() => {
        if (isMultiType) {
            if (!isArray(current_value)) throw new Error("Expected array - this should be unreachable");
            return columns.filter(c => current_value.includes(c.field));
        } else {
            return columns.find(c => c.field === current_value);
        }
    }, [current_value, columns, isMultiType]);

    return {setSelectedColumn, placeholder, type, isMultiType, columns, current_value, value};

};

/**
 * A component for selecting columns from the data store.
 * Depending on the type, this may be a single column or multiple columns, but for the purposes of this component,
 * it will only need to understand columns originating from 
 */
const ColumnDropdownComponent = observer(<T extends CTypes, M extends boolean>(gProps: ColumnSelectionProps<T, M>) => {
    const { setSelectedColumn, placeholder, isMultiType, columns, current_value, value } = useColumnDropdownValue(gProps);
    
    // should we be using a `GuiSpec<"dropdown"> | GuiSpec<"multidropdown">` here?
    // we decided against - we have some extra adornments we add and we don't want to have to deal with that in GuiSpec
    // (at least, not yet - might be tempted but there's a serious risk of types getting further out of hand)

    
    return (
        <Grid className="w-full items-center" container>
            <Grid size={"grow"}>
                <Autocomplete
                    className="w-full"
                    options={columns}
                    multiple={isMultiType}
                    //@ts-expect-error autocomplete value type seems sane - but we have an error...
                    value={value}
                    onChange={(_, value) => {
                        //! fixme <<<
                        if (!value) return; //! check if this is correct
                        if (!(isMultiType === isArray(value))) throw new Error("type mismatch");
                        if (isMultiType) {
                            // todo - need to make controlled anyway for multiple...
                            if (!isArray(value)) throw new Error("Expected array - this should be unreachable");
                            //@ts-expect-error we could use a narrower type for the setter we use here
                            setSelectedColumn(value.map(v => v.field));
                        } else {
                            if (isArray(value)) throw new Error("Unexpected array - this should be unreachable");
                            //@ts-expect-error we could use a narrower type for the setter we use here
                            setSelectedColumn(value.field);
                        }
                    }}
                    //@ts-expect-error wtf is going on, why would an individual option be an array?
                    getOptionLabel={(column) => column.name}
                    renderInput={(params) => {
                        const { key, ...p } = params as typeof params & {
                            key: string;
                        };
                        return (
                            <TextFieldExtended
                                key={key}
                                {...p}
                                placeholder={placeholder}
                            // customStartAdornment={<LinksComponent />}
                            />
                        );
                    }}
                    //@ts-expect-error wtf is going on, why would an individual option be an array? (go back to inferred type when fixed)
                    renderOption={(props, column: DataColumn<DataType>) => {
                        const { key, ...p } = props as typeof props & {
                            key: string;
                        };
                        const { datatype } = column;
                        // todo: consider an optional description prop, which we could show in a tooltip?
                        if (isMultiType) {
                            // we do hit this error... but it kinda sorta works - need to fix.
                            if (!isArray(current_value)) console.error(`Expected array - this should be unreachable (${current_value})`);
                            if (current_value === undefined) throw new Error("Expected current_value to be defined");
                            //@ts-expect-error if it was behaving, type would be narrow after throw above - not throwing because it's not behaving
                            const selected = current_value?.includes(column.field);
                            return (
                                <li key={key} {...p}>
                                    <Checkbox
                                        icon={icon}
                                        checkedIcon={checkedIcon}
                                        style={{ marginRight: 8 }}
                                        checked={selected}
                                    />
                                    {column.name}
                                    <em className="opacity-40 ml-2">({datatype})</em>
                                </li>
                            );    
                        }
                        return (
                            <li key={key} {...p}
                            onClick={(e) => {
                                p.onClick?.(e);
                            }}
                            >
                                {column.name}
                                <em className="opacity-40 ml-2">({datatype})</em>
                            </li>
                        );
                    }}
                />
            </Grid>
        </Grid>
    );
});

export default ColumnDropdownComponent;