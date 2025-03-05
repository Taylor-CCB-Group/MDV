import type React from "react";
import { useEffect, useLayoutEffect, useMemo, useState } from "react";
import { observer } from "mobx-react-lite";
import Autocomplete from "@mui/material/Autocomplete";
import { useDataStore } from "../context.js";
import type { DataColumn, DataType, FieldName, GuiSpec } from "@/charts/charts.js";
import type { Param } from "@/charts/ChartTypes.js";
import type DataStore from "@/datastore/DataStore.js";
import { g, isArray } from "@/lib/utils.js";
import { columnMatchesType } from "@/lib/columnTypeHelpers.js";
// todo - get the gui looking respectable with LinksComponent, and get it to work.
// todo - get multiple working properly.
// todo - subgroups
import LinksComponent, { RAComponent, type RowsAsColsProps } from "./LinksComponent.js";
import { TextFieldExtended } from "./TextFieldExtended.js";
import Grid from '@mui/material/Grid2';
import { Accordion, AccordionDetails, AccordionSummary, Box, Button, Paper, Tab, Tabs, Typography, useTheme } from "@mui/material";
import LinkIcon from '@mui/icons-material/Link';
import { useHighlightedForeignRows, useRowsAsColumnsLinks } from "../chartLinkHooks.js";
import type { CTypes, ColumnSelectionProps, FieldSpecs, IsMultiParam } from "@/lib/columnTypeHelpers.js";
import { inferGenericColumnGuiProps, isMultiColumn } from "@/lib/columnTypeHelpers.js";
import { action, autorun, makeAutoObservable, runInAction } from "mobx";
import { DropdownAutocompleteComponent } from "./SettingsDialogComponent.js";
import { MultiColumnQuery, RowsAsColslink, RowsAsColsQuery } from "@/links/link_utils.js";

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

    useEffect(() => {
        autorun(() => {
            if (isMultiType) {
                // todo - need to put a valid check for current value and set it accordingly
                setSelectedColumn([columns[0].field] as any);
            } else {
                if (!columns.find(c => c.field === current_value))
                    setSelectedColumn(columns[0].field as any);
            }
        });
    }, [isMultiType, columns, current_value, setSelectedColumn]);

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
const ColumnDropdown = observer(<T extends CTypes,>(gProps: ColumnSelectionProps<T> & GuiStateProps) => {
    // const props = inferGenericColumnGuiProps(gProps);
    // const { setSelectedColumn, placeholder, type } = props;
    // const { setIsAutocompleteFocused, setIsExpanded } = props;
    
    // // - this is starting to do the right thing, still massively confusing
    // const isMultiType = isMultiColProp(props);
    // const dataStore = useDataStore(props.dataStore);
    // // todo column groups
    // const columns: DataColumn<DataType>[] = useMemo(
    //     () => dataStore.columns
    //         .filter((c) => !props.exclude?.includes(c.name))
    //         .filter((c) => !c.field.includes("|")) //exclude linked columns the hacky way for now
    //         //@ts- expect-error << looks like we don't need this any more.
    //         .filter((c) => columnMatchesType(c, type))
    //         ,
    //     [dataStore, props.exclude, type],
    // );
    
    // useEffect(() => {
    //     console.log("useEffect all");
    // });

    const { setSelectedColumn, placeholder, type, setIsAutocompleteFocused, setIsExpanded,isMultiType, columns, current_value } = useColumnDropdownValue(gProps);
    
    // const linkProps = useRowsAsColumnsLinks(); //todo: arbitrary number of links, and subgroups within links
    // if (current_value && typeof current_value !== "string" && !Array.isArray(current_value)) {
    //     return <RAComponent {...props} />;
    // };
    return (
        <Grid className="w-full items-center" container>
            <Grid size={"grow"}>
                {current_value &&
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
}
            </Grid>
        </Grid>
    );
});

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
    console.log(targetColumn);
    const formattedValues = targetColumn.values.map((value, index) => (`${Object.keys(subgroups)[0]}|${value} (${Object.keys(subgroups)[0]})|${index}`))
    
    return {values: targetColumn.values, formattedValues}; //this output should be in the 'FieldName' form with pipes
}

const LinkTo = observer(<T extends CTypes,>(props: RowsAsColsProps<T>) => {
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
        type: isMultiType ? 'multidropdown' : 'dropdown',
        // name: name_column,
        label: `specific '${name}' column`, //todo different label for multiple
        values: [values],
        //this is not what we want to show in a dropdown... this is what a component will be fed if it has opted for 'active selection' mode
        current_value: [values[0]],
        func: action((v) => {
            // we don't set current_value here, we use the setSelectedColumn function
            // props.current_value = isMultiType ? v : v[0];
            // const value = `${Object.keys(subgroups)[0]}|${v[0]} (${Object.keys(subgroups)[0]})|${0}`;
            //@ts-expect-error maybe we'll make it so that we have a hook returning a setSelectedLinkedColumn function of narrower type
            setSelectedColumn(v[0])
        })
    })), [values, name, isMultiType, setSelectedColumn, props.current_value]);
    // const { current_value } = spec;
    return (
        <div className="flex flex-col" style={{ textAlign: 'left' }}>
            {/* really don't want to have this typecast here, what a nuisance! (also not really correct) */}
            <DropdownAutocompleteComponent props={spec as GuiSpec<"dropdown"> | GuiSpec<"multidropdown">} />
        </div>
    )
})

const ActiveLink = observer(<T extends CTypes,>(props: RowsAsColsProps<T>) => {
    const { linkedDs, link } = props;
    const { setSelectedColumn } = props;
    const rowNames = useHighlightedForeignRows().map(r => r.fieldName);
    const maxItems = props.multiple ? 10 : 1;

    useEffect(() => {
        //@ts-expect-error probably need isMulti logic here - may be able to refactor that into a hook
        setSelectedColumn(new RowsAsColsQuery(link, linkedDs.name, maxItems))
    }, [link, linkedDs, maxItems, setSelectedColumn]);

    return (
        <div>
            <RAComponent {...props} />
            {link.observableFields[0].fieldName}
        </div>
    );
})

/**
 * A component for selecting a column from the data store.
 * Must be in a context where `useDataStore` is available
 * (e.g. if we're in a more global dialog etc rather than a chart context,
 * this would be ambiguous).
 */
const ColumnSelectionComponent = observer(<T extends CTypes,>(props: ColumnSelectionProps<T>) => { //GuiSpec<"column">) => {
    const [isExpanded, setIsExpanded] = useState(false);
    const [isAutocompleteFocused, setIsAutocompleteFocused] = useState(false);
    const [activeTab, setActiveTab] = useState(0);
    const theme = useTheme();

    const guiProps = { isExpanded, setIsExpanded, isAutocompleteFocused, setIsAutocompleteFocused };
    const linkProps = useRowsAsColumnsLinks(); //todo: arbitrary number of links
    const rowLinkProps = useRowsAsColumnsLinks()[0]; //todo: arbitrary number of links
    //@ts-ignore nonsense about setBoolean type
    if (!linkProps) return <ColumnDropdown {...props} {...guiProps} />;
    // In general, this should (probably) be using our "(multi)dropdown" component
    // and the <LinksComponent /> can select which column options to show.

    return (
        <>
            <Paper className="mx-auto mt-8 px-4 py-2 w-full" variant="outlined">
                <div className="w-full flex justify-around text-xs font-medium">
                    <button
                        onClick={() => setActiveTab(0)}
                        type="button"
                        className="p-2 text-center border-b-2 transition-colors w-full"
                        style={{
                            borderColor: activeTab === 0 ? theme.palette.primary.main : theme.palette.divider,
                            color:
                            activeTab === 0
                                ? theme.palette.primary.main
                                : theme.palette.text.primary,
                        }}
                    >
                    Value
                    </button>
                    {rowLinkProps && (<button
                        onClick={() => setActiveTab(1)}
                        type="button"
                        className="p-2 text-center border-b-2 transition-colors w-full"
                        style={{
                        borderColor: activeTab === 1 ? theme.palette.primary.main : theme.palette.divider,
                        color:
                            activeTab === 1
                            ? theme.palette.primary.main
                            : theme.palette.text.primary,
                        }}
                    >
                    Link
                    </button>
                    )}

                    {rowLinkProps && (<button
                        onClick={() => setActiveTab(2)}
                        type="button"
                        className="p-2 text-center border-b-2 transition-colors w-full"
                        style={{
                        borderColor: activeTab === 2 ? theme.palette.primary.main : theme.palette.divider,
                        color:
                            activeTab === 2
                            ? theme.palette.primary.main
                            : theme.palette.text.primary,
                        }}
                    >
                    Active Link
                    </button>
                    )}
                </div>

                <div className="py-4">
                    {activeTab === 0 && (
                    /* we may want to show something different, especially if special value is selected... */                       
                    /* @ts-ignore setExpanded type */
                    <ColumnDropdown {...props} {...guiProps} />
                    )}
                    {activeTab === 1 && rowLinkProps && (
                            <div><LinkTo {...rowLinkProps} {...props} /></div>
                        
                    )}
                    {activeTab === 2 && rowLinkProps && (
                            <div><ActiveLink {...rowLinkProps} {...props} /></div>
                    )}
                </div>
            </Paper>
        </>
    );
    // return (
    //     <>
    //     <Accordion expanded={isExpanded} onChange={e => {
    //         if (!isAutocompleteFocused) {
    //             setIsExpanded(prev => !prev);
    //         }
    //     }}>
    //         <AccordionSummary expandIcon={<LinkIcon sx={{marginLeft: '0.2em'}} />} sx={{padding: '0 0.5em'}}>
    //             <Grid className="w-full items-center" container>
    //                 <Grid size={"grow"}>
    //                     {/* we may want to show something different, especially if special value is selected... */}                        
    //                     {/* @ts-ignore setExpanded type */}
    //                     <ColumnDropdown {...props} {...guiProps} />
    //                 </Grid>
    //             </Grid>
    //         </AccordionSummary>
    //         <AccordionDetails>
    //             <LinksComponent {...props} />
    //         </AccordionDetails>
    //     </Accordion>
    //     </>
    // );
});
export default ColumnSelectionComponent;
