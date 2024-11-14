import { observer } from "mobx-react-lite";
import { useDataStore } from "../context";
import { action, makeAutoObservable } from "mobx";
import { useHighlightedForeignRows, useRowsAsColumnsLinks } from "../chartLinkHooks";
import type { DataColumn, DataType, GuiSpec } from "@/charts/charts";
import { useMemo } from "react";
import { DropdownAutocompleteComponent } from "./SettingsDialogComponent";
import { Button, Paper } from "@mui/material";
import { isMultiColumn, type CTypes, type ColumnSelectionProps } from "@/lib/columnTypeHelpers";
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
    const isMultiType = isMultiColumn(props.type);
    // potential symbols for live link ➤ ⌁ ⇢ ⍆ ⚡︎ ► ◎ ▷ ☑︎ ⦿
    const liveSelectionName = `⦿⌁ active '${name}' selection`;
    // what do we really want here? ColumnSelectionDialogs all the way down?
    const spec: GuiSpec<'multidropdown' | 'dropdown'> = useMemo(() => makeAutoObservable({
        type: 'multidropdown',
        name: name_column,
        label: `specific '${name}' column`, //todo different label for multiple
        // I don't think we want to prepend option to dropdown - we should have a different way of showing this
        values: [targetColumn.values],
        //this is not what we want to show in a dropdown... this is what a component will be fed if it has opted for 'active selection' mode
        // current_value: props.current_value, 
        func: action((v) => {
            props.current_value = v as any;
        })
    }), [targetColumn, name_column, name]);
    const { current_value } = spec;
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

export const RAComponent = observer(<T extends CTypes,>(props: ColumnSelectionProps<T>) => {
    const { current_value } = props; //may want more props for dealing with accordion mouse event noise
    // this type check will need to change in future if we have other kinds of virtual 'query' columns
    if (!(current_value instanceof RowsAsColsQuery)) {
        console.error("Expected RowsAsColsQuery");
        return <div>Error! expected RowsAsColsQuery</div>;
    }
    const { link, linkedDsName, maxItems } = current_value;
    const multiple = isMultiColumn(props.type);
    return <Paper 
    onClick={action((e) => {
        current_value.maxItems = 100;
        e.stopPropagation();
    })}
    >⦿⌁{maxItems} '{link.name}'</Paper>;
})


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
