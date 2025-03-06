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

/**
 * A component for selecting a column from the data store - which can entail quite complex UI
 * and state for things like synthesising columns from linked data etc.
 * 
 * Must be in a context where `useDataStore` is available
 * (e.g. if we're in a more global dialog etc rather than a chart context,
 * this would be ambiguous).
 */
const ColumnSelectionComponent = observer(<T extends CTypes,>(props: ColumnSelectionProps<T>) => { //GuiSpec<"column">) => {
    const [isExpanded, setIsExpanded] = useState(false);
    const [isAutocompleteFocused, setIsAutocompleteFocused] = useState(false);
    const [activeTab, setActiveTab] = useState(0);
    const theme = useTheme();
    const { current_value } = props;

    //todo find a way to check for the type of current value and set active tab
    // useEffect(() => {
    //     // Setting the active tab based on the type of current value
    //     if (current_value) {
    //         if (typeof current_value === 'string') {
    //             if (current_value.includes("|")) setActiveTab(1);
    //             else setActiveTab(0);
    //         } else {
    //             setActiveTab(2);
    //         }
    //     }
    // }, [current_value, setActiveTab]);

    const guiProps = { isExpanded, setIsExpanded, isAutocompleteFocused, setIsAutocompleteFocused };
    const linkProps = useRowsAsColumnsLinks(); //todo: arbitrary number of links
    const rowLinkProps = useRowsAsColumnsLinks()[0]; //todo: arbitrary number of links
    //@ts-ignore nonsense about setBoolean type
    if (!linkProps) return <ColumnDropdownComponent {...props} {...guiProps} />;
    // In general, this should (probably) be using our "(multi)dropdown" component
    // and the <LinksComponent /> can select which column options to show.

    return (
        <>
        {rowLinkProps ? (
            <Paper className="mx-auto px-4 py-2 w-full" variant="outlined">
                        <div className="w-full flex justify-around text-xs font-medium">
                        <button
                            onClick={() => setActiveTab(0)}
                            type="button"
                            className="p-2 text-center border-b-2 transition-colors w-full"
                            style={{
                                borderColor: activeTab === 0 ? theme.palette.primary.main : theme.palette.divider,
                                color:
                                activeTab === 0
                                    ? theme.palette.primary.main
                                    : theme.palette.text.primary,
                            }}
                        >
                        Column
                        </button>
                        <button
                            onClick={() => setActiveTab(1)}
                            type="button"
                            className="p-2 text-center border-b-2 transition-colors w-full"
                            style={{
                            borderColor: activeTab === 1 ? theme.palette.primary.main : theme.palette.divider,
                            color:
                                activeTab === 1
                                ? theme.palette.primary.main
                                : theme.palette.text.primary,
                            }}
                        >
                        Link
                        </button>

                        <button
                            onClick={() => setActiveTab(2)}
                            type="button"
                            className="p-2 text-center border-b-2 transition-colors w-full"
                            style={{
                            borderColor: activeTab === 2 ? theme.palette.primary.main : theme.palette.divider,
                            color:
                                activeTab === 2
                                ? theme.palette.primary.main
                                : theme.palette.text.primary,
                            }}
                        >
                        Active Link
                        </button>
                    </div>

                    <div className="py-4">
                        {activeTab === 0 && (
                        /* we may want to show something different, especially if special value is selected... */                       
                        /* @ts-ignore setExpanded type */
                        <ColumnDropdownComponent {...props} {...guiProps} />
                        )}
                        {activeTab === 1 && (
                                <div><LinkToColumnComponent {...rowLinkProps} {...props} /></div>
                            
                        )}
                        {activeTab === 2 && (
                                <div><ActiveLinkComponent {...rowLinkProps} {...props} /></div>
                        )}
                    </div>
            </Paper>
        ) : (
                /* @ts-ignore setExpanded type */
                <ColumnDropdownComponent {...props} {...guiProps} />
        )}
        </>
    );
});
export default ColumnSelectionComponent;
