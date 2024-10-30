import { observer } from "mobx-react-lite";
import { useDataStore } from "../context";
import type { useDataSources } from "../hooks";
import { makeAutoObservable, observable } from "mobx";
import JsonView from "react18-json-view";
import { useRowsAsColumnsLinks } from "../chartLinkHooks";
import type { DataColumn, DataType, GuiSpec } from "@/charts/charts";
import { useMemo } from "react";
import { DropdownAutocompleteComponent } from "./SettingsDialogComponent";

type DataSource = ReturnType<typeof useDataSources>[0];
type RowsAsColsProps = ReturnType<typeof useRowsAsColumnsLinks>;

const RowsAsCols = observer(({linkedDs, rowsAsColumns} : RowsAsColsProps) => {
    const { name_column, name, subgroups } = rowsAsColumns;
    // const dataSources = useDataSources();
    const cm = window.mdv.chartManager;
    const targetColumn = cm.getDataSource(linkedDs.name).columnIndex[name_column] as DataColumn<DataType>;
    const ds = useDataStore();
    const spec: GuiSpec<'multidropdown'> = useMemo(() => makeAutoObservable({
        type: 'multidropdown',
        name: name_column,
        label: name,
        values: [targetColumn.values],
        current_value: targetColumn.values[0],
    }), [targetColumn, name_column, name]);
    const { current_value } = spec;
    const v = Array.isArray(current_value) ? current_value : [current_value];
    return (
        <div>
            <h3><em>'{linkedDs.name}'</em> rows as <em>'{ds.name}'</em> columns</h3>
            <DropdownAutocompleteComponent props={spec} />
            {v.join(', ')}
        </div>
    )
});


export default observer(function LinksComponent() {
    const link = useRowsAsColumnsLinks();
    if (!link) return null;
    const { linkedDs, rowsAsColumns } = link;
    return (
        <div>
            <JsonView src={rowsAsColumns} />
            <RowsAsCols {...link} />
        </div>
    )
});
