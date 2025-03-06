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

type setBoolean = ReturnType<typeof useState<boolean>>[1];
type GuiStateProps = {
    isExpanded: boolean;
    setIsExpanded: setBoolean;
    isAutocompleteFocused: boolean;
    setIsAutocompleteFocused: setBoolean;
}
// todo - all kinds of things to do here...
// - fix multi columns
// - column groups
// - think about how to present the links, including selecting subgroups etc.
// - maybe when more complex UI is involved, this should be displayed in a dialog.

/**
 * Check if a column selection prop is for multiple columns...
 * I do wish that there weren't so many layers of indirection etc here.
 * !Probably one of the most confusing parts of the codebase at the moment.
 */
function isMultiColProp(p: ColumnSelectionProps<any>): boolean {
    if (p.multiple === false) return false;
    return p.multiple || p.type && isMultiColumn(p.type);
}

const useColumnDropdownValue = <T extends CTypes,>(gProps: ColumnSelectionProps<T> & GuiStateProps) => {

    const props = inferGenericColumnGuiProps(gProps);
    const { setSelectedColumn, placeholder, type, current_value } = props;
    const { setIsAutocompleteFocused, setIsExpanded } = props;
    
    // - this is starting to do the right thing, still massively confusing
    //@ts-expect-error type of setSelectedColumn is wrong here
    const isMultiType = isMultiColProp(props);
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

    return {setSelectedColumn, placeholder, type, setIsAutocompleteFocused, setIsExpanded,isMultiType, columns, current_value};

};

/**
 * A component for selecting columns from the data store.
 * Depending on the type, this may be a single column or multiple columns.
 * As well as concrete columns, these columns may be 'virtual' columns representing properties that may
 * change dynamically (e.g. based on selection in a linked data source).
 * Must be in a context where `useDataStore` is available
 * (e.g. if we're in a more global dialog etc rather than a chart context,
 * this would be ambiguous).
 */
const ColumnDropdownComponent = observer(<T extends CTypes,>(gProps: ColumnSelectionProps<T> & GuiStateProps) => {
    const { setSelectedColumn, placeholder, type, setIsAutocompleteFocused, setIsExpanded,isMultiType, columns, current_value } = useColumnDropdownValue(gProps);
    
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
                                // onClick={(e) => e.stopPropagation()}
                                onFocus={() => setIsAutocompleteFocused(true)}
                                onBlur={() => setIsAutocompleteFocused(false)}
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
                                setIsExpanded(prev => prev);
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