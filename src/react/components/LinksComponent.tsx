import { observer } from "mobx-react-lite";
import { useDataStore } from "../context";
import { makeAutoObservable } from "mobx";
import { useHighlightedForeignRows, useRowsAsColumnsLinks } from "../chartLinkHooks";
import type { DataColumn, DataType, GuiSpec } from "@/charts/charts";
import { useMemo, useState } from "react";
import { DropdownAutocompleteComponent } from "./SettingsDialogComponent";
import LinkIcon from '@mui/icons-material/Link';
import { IconButton } from "@mui/material";
import { g } from "@/lib/utils";

type RowsAsColsProps = NonNullable<ReturnType<typeof useRowsAsColumnsLinks>>[0];

const RowsAsCols = observer((props : RowsAsColsProps) => {
    //! no this breaks rules of hooks
    if (!props) return null;
    const { linkedDs, link } = props;
    const [expanded, setExpanded] = useState(false);
    const rowNames = useHighlightedForeignRows().map(r => r.value);
    const { name_column, name, subgroups } = link;
    // const dataSources = useDataSources();
    const cm = window.mdv.chartManager;
    const targetColumn = cm.getDataSource(linkedDs.name).columnIndex[name_column] as DataColumn<DataType>;
    const ds = useDataStore();
    // potential symbols for live link ➤ ⌁ ⇢ ⍆ ⚡︎ ► ◎ ▷ ☑︎ ⦿
    const liveSelectionName = `⦿⌁ active '${name}' selection`;
    const spec: GuiSpec<'multidropdown'> = useMemo(() => makeAutoObservable(g({
        type: 'multidropdown',
        // name: name_column,
        label: `specific '${name}' column`, //todo different label for multiple
        // I don't think we want to prepend option to dropdown - we should have a different way of showing this
        values: [targetColumn.values],
        //this is not what we want to show in a dropdown... this is what a component will be fed if it has opted for 'active selection' mode
        current_value: rowNames, 
        func: (v) => {
            
        }
    })), [targetColumn, name_column, name, rowNames]);
    const { current_value } = spec;
    const v = Array.isArray(current_value) ? current_value : [current_value];
    return (
        <>
            <IconButton onClick={() => setExpanded(!expanded)}>
                <LinkIcon />
            </IconButton>
                {/* <h3><em>'{linkedDs.name}'</em> rows as <em>'{ds.name}'</em> columns</h3> */}
                {expanded && <DropdownAutocompleteComponent props={spec} />}

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
