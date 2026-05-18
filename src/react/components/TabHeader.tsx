import { useTheme } from "@mui/material";
import clsx from "clsx";

/** Filled circle shown before the tab label for the tab that matches `indicatorTab` (e.g. chart’s committed choice). */
const CHART_COMMIT_MARK = "\u25CF";

type TabHeaderProps<T extends string> = {
    activeTab: T;
    setActiveTab: (tab: T) => void;
    tabs: T[];
    /** When set, must match `tabs.length`; shown as button text instead of the tab id. */
    tabLabels?: string[];
    /**
     * Tab id that reflects what is already applied on the chart (e.g. committed subgroup or datasource link).
     * Shown with a leading ● so it stays readable without extra copy; may differ from `activeTab` while browsing.
     */
    indicatorTab?: T;
};

const TabHeader = <T extends string>({
    activeTab,
    setActiveTab,
    tabs,
    tabLabels,
    indicatorTab,
}: TabHeaderProps<T>) => {
    const theme = useTheme();
    // TODO consider re-arranging this so we have more proper association with tab panels etc.
    return (
        <div className="w-full flex justify-around text-xs" role="tablist">
            {tabs.map((tab, i) => {
                const selected = activeTab === tab;
                const indicatorInList =
                    indicatorTab !== undefined && tabs.includes(indicatorTab);
                const showCommitMark = indicatorInList && indicatorTab === tab;
                const labelText = tabLabels?.[i] ?? tab;
                const tabAccessibleName = showCommitMark ? `${labelText}, in active state` : labelText;
                return (
                    <button
                        key={tab}
                        onClick={() => setActiveTab(tab)}
                        type="button"
                        role="tab"
                        aria-selected={selected}
                        aria-label={tabAccessibleName}
                        className={clsx(
                            "p-2 text-center border-b-2 transition-colors w-full",
                            selected ? "font-bold" : "font-light",
                        )}
                        style={{
                            borderColor: selected ? theme.palette.primary.main : theme.palette.divider,
                            color: selected
                                ? theme.palette.primary.main
                                : theme.palette.text.primary,
                        }}
                    >
                        <span className="inline-flex items-center justify-center gap-0.5 max-w-full">
                            {showCommitMark ? (
                                <span
                                    className="shrink-0 leading-none"
                                    style={{ color: theme.palette.primary.main }}
                                    aria-hidden
                                >
                                    {CHART_COMMIT_MARK}
                                </span>
                            ) : null}
                            <span className="min-w-0 truncate">{labelText}</span>
                        </span>
                    </button>
                );
            })}
        </div>
    );
};

export default TabHeader;
