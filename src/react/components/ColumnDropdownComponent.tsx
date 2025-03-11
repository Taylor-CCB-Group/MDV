import { columnMatchesType, inferGenericColumnGuiProps } from "@/lib/columnTypeHelpers";
import type { ColumnSelectionProps, CTypes } from "@/lib/columnTypeHelpers";
import { useDataStore } from "../context";
import type { DataColumn, DataType } from "@/charts/charts";
import { useMemo } from "react";
import { observer } from "mobx-react-lite";
import { Autocomplete, Box, Checkbox, Chip } from "@mui/material";
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
            if (!isArray(current_value)) {
                // throw new Error("Expected array - this should be unreachable");
                console.error("Expected array - this should be unreachable (minor)");
                return [];
            }
            return columns.filter(c => current_value.includes(c.field));
        } else {
            return columns.find(c => c.field === current_value);
        }
    }, [current_value, columns, isMultiType]);
    console.log("value", value);

    const setValue = useMemo(() => {
        if (isMultiType) {
            return (v: DataColumn<DataType>[]) => {
                if (!isArray(v)) throw new Error("Expected array - this should be unreachable");
                //@ts-ignore kicking the can down the road, maybe a new typescript version will fix this
                setSelectedColumn(v.map(c => c.field));
            }
        } else {
            return (v: DataColumn<DataType>) => {
                if (isMultiType) throw new Error("Unexpected single column value for multi column dropdown");
                //@ts-ignore kicking the can down the road, maybe a new typescript version will fix this
                setSelectedColumn(v.field);
            }
        }
    }, [setSelectedColumn, isMultiType]);

    return {placeholder, type, isMultiType, columns, value, setValue};

};

/**
 * A component for selecting columns from the data store.
 * Depending on the type, this may be a single column or multiple columns, but for the purposes of this component,
 * it will only need to understand columns originating from 
 */
const ColumnDropdownComponent = observer(<T extends CTypes, M extends boolean>(gProps: ColumnSelectionProps<T, M>) => {
    const { placeholder, isMultiType, columns, value, setValue } = useColumnDropdownValue(gProps);
    
    // should we be using a `GuiSpec<"dropdown"> | GuiSpec<"multidropdown">` here?
    // we decided against - we have some extra adornments we add and we don't want to have to deal with that in GuiSpec
    // (at least, not yet - might be tempted but there's a serious risk of types getting further out of hand)
    type ColumnOption = DataColumn<DataType>;
    type ColumnValue = M extends true ? ColumnOption[] : ColumnOption | null;
    
    return (
        <Grid className="w-full items-center" container>
            <Grid className={"grow"}>
                <Autocomplete
                    className="w-full"
                    options={columns}
                    multiple={isMultiType}
                    value={value as ColumnValue}
                    disableCloseOnSelect
                    onChange={(_, value) => {
                        //! now that we say `value as ColumnValue`, MUI thinks the value for onChange is not an array - which it could be.
                        // *sigh*. time to move on.
                        if (!value) return; //! check if this is correct
                        if (!(isMultiType === isArray(value))) throw new Error("type mismatch");
                        //@ts-ignore kicking the can down the road, maybe a new typescript version will fix this
                        setValue(value);
                    }}
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
                            />
                        );
                    }}
                    renderTags={(value, getTagProps) => (
                        <Box
                          sx={{
                            display: "flex",
                            flexWrap: "wrap",
                            width: "100%",
                          }}
                        >
                          {value.map((option, index) => (
                            <Chip
                              label={option.name}
                              {...getTagProps({ index })}
                            />
                          ))}
                        </Box>
                    )}
                    renderOption={(props, column: DataColumn<DataType>) => {
                        const { key, ...p } = props as typeof props & {
                            key: string;
                        };
                        const { datatype } = column;
                        // todo: consider an optional description prop, which we could show in a tooltip?
                        if (isMultiType) {
                            if (!isArray(value)) throw new Error("Expected array - this should be unreachable");
                            if (value === undefined) throw new Error("Expected value to be defined");
                            const selected = value?.includes(column);
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