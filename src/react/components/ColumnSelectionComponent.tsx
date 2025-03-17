import { useEffect, useMemo, useState } from "react";
import { observer } from "mobx-react-lite";
// todo - get the gui looking respectable with LinksComponent, and get it to work.
// todo - get multiple working properly.
// todo - subgroups
import { Paper, useTheme } from "@mui/material";
import { useRowsAsColumnsLinks } from "../chartLinkHooks.js";
import type { CTypes, ColumnSelectionProps } from "@/lib/columnTypeHelpers.js";
import { paramAcceptsNumeric } from "@/lib/columnTypeHelpers";
import ColumnDropdownComponent from "./ColumnDropdownComponent.js";
import LinkToColumnComponent from "./LinkToColumnComponent.js";
import ActiveLinkComponent from "./ActiveLinkComponent.js";
import clsx from "clsx";
import { ErrorBoundary } from "react-error-boundary";
import ErrorComponentReactWrapper from "./ErrorComponentReactWrapper.js";
import { isArray } from "@/lib/utils.js";


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
 * Currently, the only type of 'link' that we are referring to is a `RowsAsCols` link.
 * That is only compatible with numeric types.
 */
function useIsLinkCompatible<T extends CTypes, M extends boolean>(props: ColumnSelectionProps<T, M>) {
    const numeric = paramAcceptsNumeric(props.type);
    const linkProps = useRowsAsColumnsLinks();
    const rowLinkProps = linkProps[0];
    return (numeric && rowLinkProps) || false;
}

function getActiveMode<T extends CTypes, M extends boolean>(props: ColumnSelectionProps<T, M>): ColumnMode {
    const { current_value } = props;
    if (current_value) {
        if (props.multiple !== isArray(current_value)) {
            throw new Error("Multiple type mismatch");
        }
        const v = isArray(current_value) ? current_value[0] : current_value;
        if (typeof v === 'string') {
            if (v.includes("|")) return "link";
            else return "column";
        } else {
            return "active link";
        }
    }
    return "column";
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
const ColumnSelectionComponent = observer(<T extends CTypes, M extends boolean>(props: ColumnSelectionProps<T, M>) => {
    const [activeTab, setActiveTab] = useState<ColumnMode>(getActiveMode(props));

    //todo check for the type of current value and set active tab
    const activeMode = useMemo(() => getActiveMode(props), [props]);
    const [initialActiveTab] = useState<ColumnMode>(activeMode);
    useEffect(() => {
        //slight annoying glitch when first mounted - useLayoutEffect didn't help.
        setActiveTab(initialActiveTab);
    }, [initialActiveTab]);


    const iIsLinkCompatible = useIsLinkCompatible(props);
    if (!iIsLinkCompatible) return <ColumnDropdownComponent {...props} />;

    return (
        <>
            {iIsLinkCompatible ? (
                <Paper className="mx-auto w-full" variant="outlined" sx={{ backgroundColor: "transparent" }}>
                    <div className="w-full flex justify-around text-xs font-light">
                        <TabHeader activeMode={activeMode} activeTab={activeTab} setActiveTab={setActiveTab} tabName="column" />
                        <TabHeader activeMode={activeMode} activeTab={activeTab} setActiveTab={setActiveTab} tabName="link" />
                        <TabHeader activeMode={activeMode} activeTab={activeTab} setActiveTab={setActiveTab} tabName="active link" />
                    </div>
                    <ErrorBoundary FallbackComponent={({error}) => <ErrorComponentReactWrapper error={error} title={error.toString()} />}>
                        <div className="p-4">
                            {activeTab === "column" && (
                                <ColumnDropdownComponent {...props} />
                            )}
                            {activeTab === "link" && (
                                <LinkToColumnComponent {...props} />
                            )}
                            {activeTab === "active link" && (
                                <ActiveLinkComponent {...props} />
                            )}
                        </div>
                    </ErrorBoundary>
                </Paper>
            ) : (
                <ColumnDropdownComponent {...props} />
            )}
        </>
    );
});
export default ColumnSelectionComponent;
