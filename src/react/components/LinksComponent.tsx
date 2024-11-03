import { observer } from "mobx-react-lite";
import { useDataStore } from "../context";
import type { useDataSources } from "../hooks";
import { makeAutoObservable, observable } from "mobx";
import JsonView from "react18-json-view";
import { useHighlightedForeignRows, useRowsAsColumnsLinks } from "../chartLinkHooks";
import type { DataColumn, DataType, GuiSpec } from "@/charts/charts";
import { useMemo, useState } from "react";
import { DropdownAutocompleteComponent } from "./SettingsDialogComponent";
import LinkIcon from '@mui/icons-material/Link';
import { Dialog, IconButton } from "@mui/material";

type RowsAsColsProps = ReturnType<typeof useRowsAsColumnsLinks>;

const RowsAsCols = observer(({linkedDs, link} : RowsAsColsProps) => {
    const [expanded, setExpanded] = useState(false);
    const rowNames = useHighlightedForeignRows().map(r => r.value);
    const { name_column, name, subgroups } = link;
    // const dataSources = useDataSources();
    const cm = window.mdv.chartManager;
    const targetColumn = cm.getDataSource(linkedDs.name).columnIndex[name_column] as DataColumn<DataType>;
    const ds = useDataStore();
    const spec: GuiSpec<'multidropdown'> = useMemo(() => makeAutoObservable({
        type: 'multidropdown',
        name: name_column,
        label: name,
        values: [["<<live link>>", ...targetColumn.values]],
        // current_value: targetColumn.values[0],
        current_value: rowNames,
    }), [targetColumn, name_column, name, rowNames]);
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
    const linkProps = useRowsAsColumnsLinks(); //todo: arbitrary number of links
    if (!linkProps) return null;
    const { linkedDs, link } = linkProps;
    return (
        <>
            {/* <JsonView src={link} /> */}
            <RowsAsCols {...linkProps} />
        </>
    )
});
