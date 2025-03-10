import { observer } from "mobx-react-lite";
import { action, makeAutoObservable, runInAction } from "mobx";
import { useHighlightedForeignRows, useRowsAsColumnsLinks } from "../chartLinkHooks";
import type { DataColumn, DataType } from "@/charts/charts";
import { useMemo } from "react";
import { DropdownAutocompleteComponent } from "./SettingsDialogComponent";
import { Button } from "@mui/material";
import { isMultiColumn, type CTypes, type ColumnSelectionProps } from "@/lib/columnTypeHelpers";
import { RowsAsColsQuery } from "@/links/link_utils";
import { g } from "@/lib/utils";

export type RowsAsColsProps<T extends CTypes, M extends boolean> = ReturnType<typeof useRowsAsColumnsLinks>[0] & ColumnSelectionProps<T, M>;


const RowsAsCols = observer(<T extends CTypes, M extends boolean>(props : RowsAsColsProps<T, M>) => {
    const { linkedDs, link } = props;
    const { setSelectedColumn } = props;
    const rowNames = useHighlightedForeignRows().map(r => r.fieldName);
    const { name_column, name, subgroups } = link;
    // const dataSources = useDataSources();
    const cm = window.mdv.chartManager;
    const targetColumn = cm.getDataSource(linkedDs.name).columnIndex[name_column] as DataColumn<DataType>;
    //@ts-expect-error need to review isMultiType logic
    const isMultiType = isMultiColumn(props.type);
    // potential symbols for live link ➤ ⌁ ⇢ ⍆ ⚡︎ ► ◎ ▷ ☑︎ ⦿
    const liveSelectionName = `⦿⌁ active '${name}' selection`;

    // what do we really want here? ColumnSelectionDialogs all the way down?
    // we need to decide how to manage non-string values in the config
    // as of now, it appears we need to explicity state the generic type which I thought was related
    // to some of the props not being of an appropriate type - but now I'm not sure.
    // ** the type of data where dealing with is the core thing we're trying to figure out how to represent from a user perspective
    // biome-ignore lint/correctness/useExhaustiveDependencies: ! this needs review !
    const spec = useMemo(() => makeAutoObservable(g<"multidropdown">({
        type: 'multidropdown',
        // name: name_column,
        label: `specific '${name}' column`, //todo different label for multiple
        // I don't think we want to prepend option to dropdown - we should have a different way of showing this
        values: [targetColumn.values],
        //this is not what we want to show in a dropdown... this is what a component will be fed if it has opted for 'active selection' mode
        current_value: [`${props.current_value}`],
        func: action((v) => {
            //@ts-expect-error need to handle multiple
            props.current_value = v[0];
        })
    })), [targetColumn, name_column, name]);
    const { current_value } = spec;
    const maxItems = props.multiple ? 10 : 1;
    return (
        <>
        {/* set the value... to something special... not just a specially formatted string */}
        <Button onClick={() => {setSelectedColumn(
            //@ts-expect-error need to handle multiple
            new RowsAsColsQuery(link, linkedDs.name, maxItems))
            }}>{liveSelectionName}</Button>
        <DropdownAutocompleteComponent props={spec} />
        </>
    )
});

export const RAComponent = observer(<T extends CTypes, M extends boolean>(props: ColumnSelectionProps<T, M>) => {
    const { current_value } = props; //may want more props for dealing with accordion mouse event noise
    // this type check will need to change in future if we have other kinds of virtual 'query' columns
    if (!(current_value instanceof RowsAsColsQuery)) {
        console.error("Expected RowsAsColsQuery");
        return <div>Error! expected RowsAsColsQuery</div>;
    }
    const { link, linkedDsName, maxItems } = current_value;
    //@ts-expect-error need to review isMultiType logic
    const multiple = isMultiColumn(props.type);
    return <div 
    onClick={(e) => {
        runInAction(() => {
            current_value.maxItems = 100;
            e.stopPropagation();
        });
    }}
    >⦿⌁{maxItems} '{link.name}'</div>;
})

export default observer(function LinksComponent<T extends CTypes, M extends boolean>(props: ColumnSelectionProps<T, M>) {
    const linkProps = useRowsAsColumnsLinks()[0]; //todo: arbitrary number of links
    if (!linkProps) return null;
    return (
        <>
            {/* <JsonView src={link} /> */}
            <RowsAsCols {...linkProps} {...props} />
        </>
    )
});
