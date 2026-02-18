import { columnMatchesType, inferGenericColumnSelectionProps } from "@/lib/columnTypeHelpers";
import type { ColumnSelectionProps, CTypes } from "@/lib/columnTypeHelpers";
import { useDataStore } from "../context";
import type { DataColumn, DataType } from "@/charts/charts";
import { useCallback, useMemo, useRef, useState } from "react";
import { observer } from "mobx-react-lite";
import { Autocomplete, Box, Button, Checkbox, Chip, Divider, Paper, type PaperProps } from "@mui/material";
import { isArray } from "@/lib/utils";
import { TextFieldExtended } from "./TextFieldExtended";
import Grid from '@mui/material/Grid2';
import CheckBoxOutlineBlankIcon from "@mui/icons-material/CheckBoxOutlineBlank";
import CheckBoxIcon from "@mui/icons-material/CheckBox";
import { useCloseOnIntersection, usePasteHandler } from "../hooks";

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

    const props = inferGenericColumnSelectionProps(gProps);
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
                setSelectedColumn(v?.field);
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
    const [selectAll, setSelectAll] = useState(isMultiType && isArray(value) ? columns.length === value?.length : false);
    const [open, setOpen] = useState(false);
    const ref = useRef<HTMLInputElement>(null);

    useCloseOnIntersection(ref, () => setOpen(false));
    // should we be using a `GuiSpec<"dropdown"> | GuiSpec<"multidropdown">` here?
    // we decided against - we have some extra adornments we add and we don't want to have to deal with that in GuiSpec
    // (at least, not yet - might be tempted but there's a serious risk of types getting further out of hand)
    type ColumnOption = DataColumn<DataType>;
    type ColumnValue = M extends true ? ColumnOption[] : ColumnOption | null;

    const handleSelectAll = useCallback(() => {
        if (!isMultiType) {
            console.error("expected multitype here");
            return;
        }

        if (selectAll) {
            // @ts-expect-error: Incompatible types
            setValue([]);
            setSelectAll(false);
        } else {
            // @ts-expect-error: Incompatible types
            setValue(columns);
            setSelectAll(true);
        }
    }, [selectAll, setValue, columns, isMultiType]);

    const handleValueChange = useCallback((newValue: DataColumn<DataType> | DataColumn<DataType>[] | null) => {
        if (!newValue) return;
        if (isMultiType) {
            (setValue as (v: DataColumn<DataType>[]) => void)(newValue as DataColumn<DataType>[]);
        } else {
            (setValue as (v: DataColumn<DataType>) => void)(newValue as DataColumn<DataType>)
        }
    }, [isMultiType, setValue]);

    // Paste handler
    const handlePaste = usePasteHandler({
        options: columns,
        multiple: isMultiType,
        currentValue: value || null,
        setValue: handleValueChange,
        getLabel: (column: DataColumn<DataType>) => column.name,
        getValue: (column: DataColumn<DataType>) => column,
    });
    
    return (
        <Grid className="w-full items-center" container>
            <Grid className={"grow"}>
                <Autocomplete
                    //! A future use-case can be that the columns will be added at runtime, but MUI's autocomplete has an existing issue. Refer: https://github.com/mui/material-ui/issues/29508 
                    className="w-full"
                    options={columns}
                    multiple={isMultiType}
                    value={value as ColumnValue}
                    disableCloseOnSelect={isMultiType}
                    open={open}
                    onOpen={() => setOpen(true)}
                    onClose={() => setOpen(false)}
                    ref={ref}
                    onChange={(_, value) => {
                        //! now that we say `value as ColumnValue`, MUI thinks the value for onChange is not an array - which it could be.
                        // *sigh*. time to move on.
                        if (!value) return; //! check if this is correct
                        if (!(isMultiType === isArray(value))) throw new Error("type mismatch");
                        if (isMultiType && isArray(value) && value.length === 0) {
                            setSelectAll(false);
                        }
                        //@ts-ignore kicking the can down the road, maybe a new typescript version will fix this
                        setValue(value);
                    }}
                    getOptionLabel={(column) => column.name}
                    renderInput={(params) => {
                        const { key, InputProps, ...p } = params as typeof params & {
                            key: string;
                        };
                        return (
                            <TextFieldExtended
                                key={key}
                                {...p}
                                placeholder={placeholder}
                                slotProps={{
                                    input: {
                                        ...InputProps,
                                        
                                        onPaste: handlePaste,
                                    }
                                }}
                            />
                        );
                        
                    }}
                    renderTags={(value, getTagProps) => {
                        // custom logic to limit the tags, this is required because we are overriding the way tags are rendered
                        const limit = 5;
                        const displayed = value.slice(0, limit);
                        const remaining = value.length - limit;

                        return (
                            <Box
                                sx={{
                                    display: "flex",
                                    flexWrap: "wrap",
                                    width: "100%",
                                }}
                            >
                                {displayed.map((option, index) => {
                                    const { key, ...tagProps } = getTagProps({ index });
                                    return (
                                    <Chip
                                        key={key}
                                        label={option.name}
                                        {...tagProps}
                                    />
                                )})}
                                {remaining > 0 && <Chip label={`+${remaining}`} />}
                            </Box>
                        );
                    }}
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
                    // biome-ignore lint/correctness/useExhaustiveDependencies: not including handleSelectAll because it rerenders the paper component
                    PaperComponent={useMemo(() => (paperProps: PaperProps) => {
                        const { children, ...restPaperProps } = paperProps;
                        return (
                          <Paper {...restPaperProps}>
                            {isMultiType &&
                            (<>
                                <Box
                                    onMouseDown={(e) => e.preventDefault()}
                                    py={1}
                                    px={2}
                                >
                                    <Button onClick={handleSelectAll} color="inherit" fullWidth>
                                        {selectAll ? "Unselect All" : "Select All"}
                                    </Button>
                                </Box>
                                <Divider />
                            </>)}
                            {children}
                          </Paper>
                        );
                      }, [isMultiType, selectAll])}
                />
            </Grid>
        </Grid>
    );
});

export default ColumnDropdownComponent;