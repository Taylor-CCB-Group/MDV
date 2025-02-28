import React, { useEffect, useMemo, useState } from "react";
import { observer } from "mobx-react-lite";
import Autocomplete from "@mui/material/Autocomplete";
import { useDataStore } from "../context.js";
import type { DataColumn, DataType, FieldName } from "@/charts/charts.js";
import type { Param } from "@/charts/ChartTypes.js";
import type DataStore from "@/datastore/DataStore.js";
import { g, isArray } from "@/lib/utils.js";
import { columnMatchesType } from "@/lib/columnTypeHelpers.js";
// todo - get the gui looking respectable with LinksComponent, and get it to work.
// todo - get multiple working properly.
// todo - subgroups
import LinksComponent, { RAComponent, RowsAsColsProps } from "./LinksComponent.js";
import { TextFieldExtended } from "./TextFieldExtended.js";
import Grid from '@mui/material/Grid2';
import { Accordion, AccordionDetails, AccordionSummary, Box, Button, Paper, Tab, Tabs, Typography, useTheme } from "@mui/material";
import LinkIcon from '@mui/icons-material/Link';
import { useHighlightedForeignRows, useRowsAsColumnsLinks } from "../chartLinkHooks.js";
import type { CTypes, ColumnSelectionProps } from "@/lib/columnTypeHelpers.js";
import { inferGenericColumnGuiProps, isMultiColumn } from "@/lib/columnTypeHelpers.js";
import { action, makeAutoObservable } from "mobx";
import { DropdownAutocompleteComponent } from "./SettingsDialogComponent.js";
import { RowsAsColsQuery } from "@/links/link_utils.js";

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
    const props = inferGenericColumnGuiProps(gProps);
    const { setSelectedColumn, placeholder, type } = props;
    const { setIsAutocompleteFocused, setIsExpanded, current_value } = props;
    
    // - this is starting to do the right thing, still massively confusing
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
    // const linkProps = useRowsAsColumnsLinks(); //todo: arbitrary number of links, and subgroups within links
    // if (current_value && typeof current_value !== "string" && !Array.isArray(current_value)) {
    //     return <RAComponent {...props} />;
    // };
    return (
        <Grid className="w-full items-center" container>
            <Grid size={"grow"}>
                <Autocomplete
                    className="w-full"
                    options={columns}
                    multiple={isMultiType}
                    value={current_value ? columns.find(c => c.name === current_value) : null}
                    onChange={(_, value) => {
                        //! fixme
                        if (!value) return; //! check if this is correct
                        if (!(isMultiType === isArray(value))) throw "type mismatch";
                        if (isMultiType) {
                            // todo - need to make controlled anyway for multiple...
                            if (!isArray(value)) throw "Expected array - and really didn't expect to get here";
                            // haven't quite managed to infer the setSelectedColumn type here...
                            setSelectedColumn(value.map(v => v.field) as any);
                        } else {
                            if (isArray(value)) throw "Expected single value - and really didn't expect to get here";
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

type CustomTabPanelProps = {
    children?: React.ReactNode,
    index: number,
    value: number,
}

const CustomTabPanel = (props: CustomTabPanelProps) => {
    const {children, index, value, ...other} = props;
    return (
        <div
            role="tabpanel"
            hidden={value !== index}
            id={`tab-${index}`}
            aria-labelledby={`tab-${index}`}
            style={{padding: '1em 0.5em', width: "100%"}}
            {...other}
            >
            {value === index && children}
        </div>
    )
};

const LinkTo = observer(<T extends CTypes,>(props: RowsAsColsProps<T>) => {
    const { linkedDs, link } = props;
    const { name_column, name, subgroups } = link;
    // const dataSources = useDataSources();
    const cm = window.mdv.chartManager;
    const targetColumn = cm.getDataSource(linkedDs.name).columnIndex[name_column] as DataColumn<DataType>;
    //@ts-expect-error need to review isMultiType logic
    const isMultiType = isMultiColumn(props.type);
    // potential symbols for live link ➤ ⌁ ⇢ ⍆ ⚡︎ ► ◎ ▷ ☑︎ ⦿
    const liveSelectionName = `⦿⌁ active '${name}' selection`;
    const spec = useMemo(() => makeAutoObservable(g<"multidropdown">({
        type: 'multidropdown',
        // name: name_column,
        label: `specific '${name}' column`, //todo different label for multiple
        // I don't think we want to prepend option to dropdown - we should have a different way of showing this
        values: [targetColumn.values],
        //this is not what we want to show in a dropdown... this is what a component will be fed if it has opted for 'active selection' mode
        current_value: [`${props.current_value}`],
        func: action((v) => {
            props.current_value = v[0];
        })
    })), [targetColumn, name_column, name]);
    const { current_value } = spec;
    return (
        <DropdownAutocompleteComponent props={spec} />
    )
})

const ActiveLink = observer(<T extends CTypes,>(props: RowsAsColsProps<T>) => {
    const { linkedDs, link } = props;
    const { setSelectedColumn } = props;
    const rowNames = useHighlightedForeignRows().map(r => r.fieldName);
    const maxItems = props.multiple ? 10 : 1;

    useEffect(() => {
        setSelectedColumn(new RowsAsColsQuery(link, linkedDs.name, maxItems))
    }, [link, linkedDs, maxItems]);

    return (
        <div>
            <RAComponent {...props} />
            Selected Active Link
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
    console.log("props", props);
    const theme = useTheme();
    function a11yProps(index: number) {
        return {
          id: `tab-${index}`,
          'aria-controls': `tab-${index}`,
        };
      }

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
                <div className="border-b border-gray-200 w-full flex justify-around">
                    <button
                        onClick={() => setActiveTab(0)}
                        className="p-2 text-center border-b-2 transition-colors w-full"
                        style={{
                            borderColor: activeTab === 0 ? theme.palette.primary.main : "transparent",
                            color:
                            activeTab === 0
                                ? theme.palette.primary.main
                                : theme.palette.text.primary,
                        }}
                    >
                    Value
                    </button>
                    <button
                        onClick={() => setActiveTab(1)}
                        className="p-2 text-center border-b-2 transition-colors w-full"
                        style={{
                        borderColor: activeTab === 1 ? theme.palette.primary.main : "transparent",
                        color:
                            activeTab === 1
                            ? theme.palette.primary.main
                            : theme.palette.text.primary,
                        }}
                    >
                    Link to
                    </button>
                    <button
                        onClick={() => setActiveTab(2)}
                        className="p-2 text-center border-b-2 transition-colors w-full"
                        style={{
                        borderColor: activeTab === 2 ? theme.palette.primary.main : "transparent",
                        color:
                            activeTab === 2
                            ? theme.palette.primary.main
                            : theme.palette.text.primary,
                        }}
                    >
                    Active Link to
                    </button>
                </div>

                <div className="py-4">
                    {activeTab === 0 && (
                    /* we may want to show something different, especially if special value is selected... */                       
                    /* @ts-ignore setExpanded type */
                    <ColumnDropdown {...props} {...guiProps} />
                    )}
                    {activeTab === 1 && (
                    rowLinkProps && (
                            <div><LinkTo {...rowLinkProps} {...props} /></div>
                        )
                    )}
                    {activeTab === 2 && (
                    rowLinkProps && (
                            <div><ActiveLink {...rowLinkProps} {...props} /></div>
                        )
                    )}
                </div>
            </Paper>
        </>
        // <>
        // <Paper sx={{width: "100%"}}>
        //     <Tabs value={tabValue} onChange={(_e, value) => setTabValue(value)} variant="scrollable">
        //         <Tab value={0} label="Value" {...a11yProps} />
        //         <Tab value={1} label="Link to" {...a11yProps} />
        //         <Tab value={2} label="Active Link to" {...a11yProps} />
        //     </Tabs>
        //     <CustomTabPanel index={0} value={tabValue}>
        //         {/* we may want to show something different, especially if special value is selected... */}                        
        //         {/* @ts-ignore setExpanded type */}
        //         <ColumnDropdown {...props} {...guiProps} />
        //     </CustomTabPanel>
        //     {rowLinkProps && (
        //         <CustomTabPanel index={1} value={tabValue}>
        //             <div><LinkTo {...rowLinkProps} {...props} /></div>
        //         </CustomTabPanel>
        //     )}
        //     {rowLinkProps && (
        //         <CustomTabPanel index={2} value={tabValue}>
        //             <div><ActiveLink {...rowLinkProps} {...props} /></div>
        //         </CustomTabPanel>
        //     )}
        // </Paper>
        // </>
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
