import { columnMatchesType, ColumnSelectionProps, CTypes, inferGenericColumnGuiProps, isMultiColumn } from "@/lib/columnTypeHelpers";
import { useDataStore } from "../context";
import { DataColumn, DataType } from "@/charts/charts";
import { useEffect, useMemo, useState } from "react";
import { autorun } from "mobx";
import { observer } from "mobx-react-lite";
import { Autocomplete } from "@mui/material";
import { isArray } from "@/lib/utils";
import { TextFieldExtended } from "./TextFieldExtended";
import Grid from '@mui/material/Grid2';

// type setBoolean = ReturnType<typeof useState<boolean>>[1];
// type GuiStateProps = {
//     isExpanded: boolean;
//     setIsExpanded: setBoolean;
//     isAutocompleteFocused: boolean;
//     setIsAutocompleteFocused: setBoolean;
// }

// todo - all kinds of things to do here...
// - make sure only compatible columns are shown
// - fix multi columns
// - column groups

/**
 * Check if a column selection prop is for multiple columns...
 * I do wish that there weren't so many layers of indirection etc here.
 * !Probably one of the most confusing parts of the codebase at the moment.
 */
function isMultiColProp(p: ColumnSelectionProps<any>): boolean {
    if (p.multiple === false) return false;
    return p.multiple || p.type && isMultiColumn(p.type);
}

const useColumnDropdownValue = <T extends CTypes,>(gProps: ColumnSelectionProps<T>) => {

    const props = inferGenericColumnGuiProps(gProps);
    const { setSelectedColumn, placeholder, type, current_value } = props;
    
    // - this is starting to do the right thing, still massively confusing
    //@ts-expect-error type of setSelectedColumn is wrong here
    const isMultiType = isMultiColProp(props);
    // are we sure that props.dataStore will be right 
    // - are there any inconsistencies between this and the current DataStoreContext?
    const dataStore = useDataStore(props.dataStore);
    // todo column groups
    const columns: DataColumn<DataType>[] = useMemo(
        () => dataStore.columns
            .filter((c) => !props.exclude?.includes(c.name))
            .filter((c) => !c.field.includes("|")) //exclude linked columns the hacky way for now
            //@ts- expect-error << looks like we don't need this any more.
            .filter((c) => columnMatchesType(c, type))
            ,
        [dataStore, props.exclude, type],
    );

    // todo make a setLocalColumn function, with single / multi type, but only taking strings
    // also allow "None" if optional
    // check why SelectionDialog is not working (multi-column)

    return {setSelectedColumn, placeholder, type, isMultiType, columns, current_value};

};

/**
 * A component for selecting columns from the data store.
 * Depending on the type, this may be a single column or multiple columns, but for the purposes of this component,
 * it will only need to understand columns originating from 
 */
const ColumnDropdownComponent = observer(<T extends CTypes,>(gProps: ColumnSelectionProps<T>) => {
    const { setSelectedColumn, placeholder, type, isMultiType, columns, current_value } = useColumnDropdownValue(gProps);
    
    return (
        <Grid className="w-full items-center" container>
            <Grid size={"grow"}>
                <Autocomplete
                    className="w-full"
                    options={columns}
                    multiple={isMultiType}
                    value={columns.find(c => c.field === current_value) || null}
                    onChange={(_, value) => {
                        //! fixme
                        if (!value) return; //! check if this is correct
                        if (!(isMultiType === isArray(value))) throw new Error("type mismatch");
                        if (isMultiType) {
                            // todo - need to make controlled anyway for multiple...
                            if (!isArray(value)) throw new Error("Expected array - and really didn't expect to get here");
                            //@ts-expect-error haven't quite managed to infer the setSelectedColumn type here...
                            setSelectedColumn(value.map(v => v.field));
                        } else {
                            if (isArray(value)) throw new Error("Expected single value - and really didn't expect to get here");
                            //@ts-expect-error haven't quite managed to infer the setSelectedColumn type here...
                            setSelectedColumn(value.field);
                        }
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
                            // customStartAdornment={<LinksComponent />}
                            />
                        );
                    }}
                    renderOption={(props, column) => {
                        const { key, ...p } = props as typeof props & {
                            key: string;
                        };
                        const { datatype } = column;
                        // todo: consider an optional description prop, which we could show in a tooltip?
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