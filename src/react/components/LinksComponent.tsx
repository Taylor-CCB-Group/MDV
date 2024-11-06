import { observer } from "mobx-react-lite";
import { useDataStore } from "../context";
import { makeAutoObservable } from "mobx";
import { useHighlightedForeignRows, useRowsAsColumnsLinks } from "../chartLinkHooks";
import type { DataColumn, DataType, GuiSpec } from "@/charts/charts";
import { useMemo } from "react";
import { DropdownAutocompleteComponent } from "./SettingsDialogComponent";
import { Button, FormControl, FormControlLabel, Radio, RadioGroup } from "@mui/material";

type RowsAsColsProps = ReturnType<typeof useRowsAsColumnsLinks>[0];

const RowsAsCols = observer((props : RowsAsColsProps) => {
    const { linkedDs, link } = props;
    const rowNames = useHighlightedForeignRows().map(r => r.value);
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
        // current_value: targetColumn.values[0],
        current_value: rowNames,
    }), [targetColumn, name_column, name, rowNames]);
    const { current_value } = spec;
    const v = Array.isArray(current_value) ? current_value : [current_value];
    return (
        <>
        <Button onClick={() => {}}>{liveSelectionName}</Button>
        <DropdownAutocompleteComponent props={spec} />
        {/* <FormControl>
            <RadioGroup>
                <FormControlLabel control={<Radio size="small" />} label={liveSelectionName} />
                
            </RadioGroup>
        </FormControl> */}
        </>
    )
});


export default observer(function LinksComponent() {
    const linkProps = useRowsAsColumnsLinks()[0]; //todo: arbitrary number of links
    if (!linkProps) return null;
    return (
        <>
            {/* <JsonView src={link} /> */}
            <RowsAsCols {...linkProps} />
        </>
    )
});
