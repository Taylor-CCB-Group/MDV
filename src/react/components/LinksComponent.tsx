import { observer } from "mobx-react-lite";
import { useDataStore } from "../context";
import { makeAutoObservable } from "mobx";
import { useHighlightedForeignRows, useRowsAsColumnsLinks } from "../chartLinkHooks";
import type { DataColumn, DataType, FieldName, GuiSpec } from "@/charts/charts";
import { useMemo } from "react";
import { DropdownAutocompleteComponent } from "./SettingsDialogComponent";
import { Button, FormControl, FormControlLabel, Radio, RadioGroup } from "@mui/material";
import type { CTypes, ColumnSelectionProps } from "@/lib/columnTypeHelpers";
import { RowsAsColsQuery } from "@/links/link_utils";

type RowsAsColsProps<T extends CTypes> = ReturnType<typeof useRowsAsColumnsLinks>[0] & ColumnSelectionProps<T>;

const RowsAsCols = observer(<T extends CTypes,>(props : RowsAsColsProps<T>) => {
    const { linkedDs, link } = props;
    const { setSelectedColumn } = props;
    const rowNames = useHighlightedForeignRows().map(r => r.fieldName);
    const { name_column, name, subgroups } = link;
    // const dataSources = useDataSources();
    const cm = window.mdv.chartManager;
    const targetColumn = cm.getDataSource(linkedDs.name).columnIndex[name_column] as DataColumn<DataType>;
    const ds = useDataStore();
    // potential symbols for live link ➤ ⌁ ⇢ ⍆ ⚡︎ ► ◎ ▷ ☑︎ ⦿
    const liveSelectionName = `⦿⌁ active '${name}' selection`;
    const spec: GuiSpec<'multidropdown'> = useMemo(() => makeAutoObservable({
        type: 'multidropdown',
        name: name_column,
        label: `specific '${name}' column`, //todo different label for multiple
        // I don't think we want to prepend option to dropdown - we should have a different way of showing this
        values: [targetColumn.values],
        //this is not what we want to show in a dropdown... this is what a component will be fed if it has opted for 'active selection' mode
        current_value: [], 
        func: (v) => {
            
        }
    }), [targetColumn, name_column, name]);
    const { current_value } = spec;
    const v = Array.isArray(current_value) ? current_value : [current_value];
    const maxItems = props.multiple ? 10 : 1;
    return (
        <>
        {/* set the value... to something special... not just a specially formatted string */}
        <Button onClick={() => {setSelectedColumn(new RowsAsColsQuery(link, linkedDs.name, maxItems))}}>{liveSelectionName}</Button>
        <DropdownAutocompleteComponent props={spec} />
        {/* <FormControl>
            <RadioGroup>
                <FormControlLabel control={<Radio size="small" />} label={liveSelectionName} />
                
            </RadioGroup>
        </FormControl> */}
        </>
    )
});


export default observer(function LinksComponent<T extends CTypes,>(props: ColumnSelectionProps<T>) {
    const linkProps = useRowsAsColumnsLinks()[0]; //todo: arbitrary number of links
    if (!linkProps) return null;
    return (
        <>
            {/* <JsonView src={link} /> */}
            <RowsAsCols {...linkProps} {...props} />
        </>
    )
});
