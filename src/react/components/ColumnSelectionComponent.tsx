import { useEffect, useState } from "react";
import { observer } from "mobx-react-lite";
// todo - get the gui looking respectable with LinksComponent, and get it to work.
// todo - get multiple working properly.
// todo - subgroups
import { Paper, useTheme } from "@mui/material";
import { useRowsAsColumnsLinks } from "../chartLinkHooks.js";
import type { CTypes, ColumnSelectionProps } from "@/lib/columnTypeHelpers.js";
import ColumnDropdownComponent from "./ColumnDropdownComponent.js";
import LinkToColumnComponent from "./LinkToColumnComponent.js";
import ActiveLinkComponent from "./ActiveLinkComponent.js";


type ColumnMode = "column" | "link" | "activeLink";

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
    const theme = useTheme();
    
    //todo check for the type of current value and set active tab
    // const { current_value } = props;
    // useEffect(() => {
    //     // Setting the active tab based on the type of current value
    //     if (current_value) {
    //         if (typeof current_value === 'string') {
    //             if (current_value.includes("|")) setActiveTab("link");
    //             else setActiveTab("column");
    //         } else {
    //             setActiveTab("activeLink");
    //         }
    //     }
    // }, [current_value, setActiveTab]);

    const linkProps = useRowsAsColumnsLinks(); //todo: arbitrary number of links
    const rowLinkProps = useRowsAsColumnsLinks()[0]; //todo: arbitrary number of links
    if (!linkProps) return <ColumnDropdownComponent {...props} />;

    return (
        <>
        {rowLinkProps ? (
            <Paper className="mx-auto w-full" variant="outlined" sx={{backgroundColor: "transparent"}}>
                        <div className="w-full flex justify-around text-xs font-light">
                        <button
                            onClick={() => setActiveTab("column")}
                            type="button"
                            className="p-2 text-center border-b-2 transition-colors w-full"
                            style={{
                                borderColor: activeTab === "column" ? theme.palette.primary.main : theme.palette.divider,
                                color:
                                    activeTab === "column"
                                        ? theme.palette.primary.main
                                        : theme.palette.text.primary,
                            }}
                        >
                            Column
                        </button>
                        <button
                            onClick={() => setActiveTab("link")}
                            type="button"
                            className="p-2 text-center border-b-2 transition-colors w-full"
                            style={{
                                borderColor: activeTab === "link" ? theme.palette.primary.main : theme.palette.divider,
                                color:
                                    activeTab === "link"
                                        ? theme.palette.primary.main
                                        : theme.palette.text.primary,
                            }}
                        >
                            Link
                        </button>

                        <button
                            onClick={() => setActiveTab("activeLink")}
                            type="button"
                            className="p-2 text-center border-b-2 transition-colors w-full"
                            style={{
                                borderColor: activeTab === "activeLink" ? theme.palette.primary.main : theme.palette.divider,
                                color:
                                    activeTab === "activeLink"
                                        ? theme.palette.primary.main
                                        : theme.palette.text.primary,
                            }}
                        >
                            Active Link
                        </button>
                    </div>

                    <div className="p-4">
                        {activeTab === "column" && (
                            <ColumnDropdownComponent {...props} />
                        )}
                        {activeTab === "link" && (
                            <div><LinkToColumnComponent {...rowLinkProps} {...props} /></div>

                        )}
                        {activeTab === "activeLink" && (
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
