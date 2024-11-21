import { useMemo, useState } from "react";
import { observer } from "mobx-react-lite";
import Autocomplete from "@mui/material/Autocomplete";
import { useDataStore } from "../context.js";
import type { DataColumn, DataType, FieldName } from "@/charts/charts.js";
import type { Param } from "@/charts/ChartTypes.js";
import type DataStore from "@/datastore/DataStore.js";
import { columnMatchesType, isArray } from "@/lib/utils.js";
// todo - get the gui looking respectable with LinksComponent, and get it to work.
// todo - get multiple working properly.
// todo - subgroups
import LinksComponent, { RAComponent } from "./LinksComponent.js";
import { TextFieldExtended } from "./TextFieldExtended.js";
import Grid from '@mui/material/Grid2';
import { Accordion, AccordionDetails, AccordionSummary, Typography } from "@mui/material";
import LinkIcon from '@mui/icons-material/Link';
import { useRowsAsColumnsLinks } from "../chartLinkHooks.js";
import type { CTypes, ColumnSelectionProps } from "@/lib/columnTypeHelpers.js";
import { inferGenericColumnGuiProps, isMultiColumn } from "@/lib/columnTypeHelpers.js";

type setBoolean = ReturnType<typeof useState<boolean>>[1];
type GuiStateProps = {
    isExpanded: boolean;
    setIsExpanded: setBoolean;
    isAutocompleteFocused: boolean;
    setIsAutocompleteFocused: setBoolean;
}


/**
 * A component for selecting a column from the data store.
 * Must be in a context where `useDataStore` is available
 * (e.g. if we're in a more global dialog etc rather than a chart context,
 * this would be ambiguous).
 */
const ColumnDropdown = observer(<T extends CTypes,>(gProps: ColumnSelectionProps<T> & GuiStateProps) => {
    const props = inferGenericColumnGuiProps(gProps);
    const { setSelectedColumn, placeholder, type } = props;
    const { setIsAutocompleteFocused, setIsExpanded, current_value } = props;
    const isMultiType = isMultiColumn(type);
    // const multiple = props.multiple || isMultiType;
    const dataStore = useDataStore(props.dataStore);
    // todo column groups
    const columns: DataColumn<DataType>[] = useMemo(
        () => dataStore.columns
            .filter((c) => !props.exclude?.includes(c.name))
            //@ts-expect-error
            .filter((c) => columnMatchesType(c, type))
            ,
        [dataStore, props.exclude, type],
    );
    // const [isExpanded, setIsExpanded] = useState(false);
    // const [isAutocompleteFocused, setIsAutocompleteFocused] = useState(false);

    // const guiProps = { isExpanded, setIsExpanded, isAutocompleteFocused, setIsAutocompleteFocused };
    // const linkProps = useRowsAsColumnsLinks(); //todo: arbitrary number of links
    // if (!linkProps) return <ColumnDropdown {...props} {...guiProps} />;
    // In general, this should (probably) be using our "(multi)dropdown" component
    // and the <LinksComponent /> can select which column options to show.
    if (current_value && typeof current_value !== "string" && !Array.isArray(current_value)) {
        return <RAComponent {...props} />;
    };
    return (
        <Grid className="w-full items-center" container>
            <Grid size={"grow"}>
                <Autocomplete
                    className="w-full"
                    options={columns}
                    multiple={isMultiType}
                    // value={current_value ? columns.find(c => c.name === current_value) : null}
                    onChange={(_, value) => {
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
                                onClick={(e) => e.stopPropagation()}
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

/**
 * A component for selecting a column from the data store.
 * Must be in a context where `useDataStore` is available
 * (e.g. if we're in a more global dialog etc rather than a chart context,
 * this would be ambiguous).
 */
const ColumnSelectionComponent = observer(<T extends CTypes,>(props: ColumnSelectionProps<T>) => { //GuiSpec<"column">) => {
    const [isExpanded, setIsExpanded] = useState(false);
    const [isAutocompleteFocused, setIsAutocompleteFocused] = useState(false);

    const guiProps = { isExpanded, setIsExpanded, isAutocompleteFocused, setIsAutocompleteFocused };
    const linkProps = useRowsAsColumnsLinks(); //todo: arbitrary number of links
    //@ts-expect-error nonsense about setBoolean type
    if (!linkProps) return <ColumnDropdown {...props} {...guiProps} />;
    // In general, this should (probably) be using our "(multi)dropdown" component
    // and the <LinksComponent /> can select which column options to show.
    return (
        <>
        <Accordion expanded={isExpanded} onChange={e => {
            if (!isAutocompleteFocused) {
                setIsExpanded(prev => !prev);
            }
        }}>
            <AccordionSummary expandIcon={<LinkIcon sx={{marginLeft: '0.2em'}} />} sx={{padding: '0 0.5em'}}>
                <Grid className="w-full items-center" container>
                    <Grid size={"grow"}>
                        {/* we may want to show something different, especially if special value is selected... */}                        
                        {/* @ts-expect-error */}
                        <ColumnDropdown {...props} {...guiProps} />
                    </Grid>
                </Grid>
            </AccordionSummary>
            <AccordionDetails>
                <LinksComponent {...props} />
            </AccordionDetails>
        </Accordion>
        </>
    );
});
export default ColumnSelectionComponent;
