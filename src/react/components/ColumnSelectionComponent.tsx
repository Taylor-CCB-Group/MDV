import { useEffect, useMemo, useState } from "react";
import { observer } from "mobx-react-lite";
// todo - get the gui looking respectable with LinksComponent, and get it to work.
// todo - get multiple working properly.
// todo - subgroups
import { Paper, useTheme } from "@mui/material";
import { useRowsAsColumnsLinks } from "../chartLinkHooks.js";
import type { CTypes, ColumnSelectionProps } from "@/lib/columnTypeHelpers.js";
import { isMultiColumn } from "@/lib/columnTypeHelpers";
import ColumnDropdownComponent from "./ColumnDropdownComponent.js";
import LinkToColumnComponent from "./LinkToColumnComponent.js";
import ActiveLinkComponent from "./ActiveLinkComponent.js";
import clsx from "clsx";


type ColumnMode = "column" | "link" | "active link";

type TabHeaderProps = {
    activeTab: ColumnMode;
    setActiveTab: (tab: ColumnMode) => void;
    tabName: ColumnMode;
    activeMode?: ColumnMode;
}

const TabHeader = ({ activeTab, setActiveTab, tabName, activeMode }: TabHeaderProps) => {
    const theme = useTheme();
    const isActiveMode = activeMode === tabName;
    return (
        <button
            onClick={() => setActiveTab(tabName)}
            type="button"
            className={clsx("p-2 text-center border-b-2 transition-colors w-full", isActiveMode ? "font-bold" : "font-light")}
            style={{
                borderColor: activeTab === tabName ? theme.palette.primary.main : theme.palette.divider,
                color:
                    activeTab === tabName
                        ? theme.palette.primary.main
                        : theme.palette.text.primary,
            }}
        >
            {tabName}
        </button>
    );
}

/**
 * A component for selecting a column from the data store - which can entail quite complex UI
 * and state for things like synthesising columns from linked data etc.
 * 
 * When 'linked' data is available, an interface is shown allowing for switching between different 
 * modes of column selection, each of which has it's own UI as well as hooks for adapting the complex generic
 * props and narrowing them to simpler types related specifically to the corresponding mode.
 * 
 * As of this writing, this more complex functionality is specifically related to `RowsAsCols` links.
 * 
 * Must either be used in a context where a `DataStoreContext` is available, or provided with a `dataStore` prop.
 * (e.g. if we're in a more global dialog etc rather than a chart context, this would be ambiguous).
 */
const ColumnSelectionComponent = observer(<T extends CTypes,>(props: ColumnSelectionProps<T>) => {
    const [activeTab, setActiveTab] = useState<ColumnMode>("column");

    //todo check for the type of current value and set active tab
    const { current_value } = props;
    
    //@ts-expect-error need to review isMultiType logic
    const isMultiType = isMultiColumn(props.type);

    const activeMode = useMemo(() => {
        if (current_value) {
            //@ts-expect-error need to review isMultiType logic
            const v = isMultiType ? current_value : current_value[0];
            if (typeof v === 'string') {
                if (v.includes("|")) return "link";
                else return "column";
            } else {
                return "active link";
            }
        }
        return "column";
    }, [current_value, isMultiType]);
    const [initialActiveTab] = useState<ColumnMode>(activeMode);
    useEffect(() => {
        //slight annoying glitch when first mounted - useLayoutEffect didn't help.
        setActiveTab(initialActiveTab);
    }, [initialActiveTab]);


    const linkProps = useRowsAsColumnsLinks(); //todo: arbitrary number of links
    const rowLinkProps = useRowsAsColumnsLinks()[0]; //todo: arbitrary number of links
    if (!linkProps) return <ColumnDropdownComponent {...props} />;

    return (
        <>
            {rowLinkProps ? (
                <Paper className="mx-auto w-full" variant="outlined" sx={{ backgroundColor: "transparent" }}>
                    <div className="w-full flex justify-around text-xs font-light">
                        <TabHeader activeMode={activeMode} activeTab={activeTab} setActiveTab={setActiveTab} tabName="column" />
                        <TabHeader activeMode={activeMode} activeTab={activeTab} setActiveTab={setActiveTab} tabName="link" />
                        <TabHeader activeMode={activeMode} activeTab={activeTab} setActiveTab={setActiveTab} tabName="active link" />
                    </div>

                    <div className="p-4">
                        {activeTab === "column" && (
                            <ColumnDropdownComponent {...props} />
                        )}
                        {activeTab === "link" && (
                            <div><LinkToColumnComponent {...rowLinkProps} {...props} /></div>

                        )}
                        {activeTab === "active link" && (
                            <div><ActiveLinkComponent {...rowLinkProps} {...props} /></div>
                        )}
                    </div>
                </Paper>
            ) : (
                <ColumnDropdownComponent {...props} />
            )}
        </>
    );
});
export default ColumnSelectionComponent;
