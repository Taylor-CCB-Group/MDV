import { useTheme } from "@mui/material";
import clsx from "clsx";

type TabHeaderProps<T extends string> = {
    activeTab: T;
    setActiveTab: (tab: T) => void;
    tabs: T[];
};

const TabHeader = <T extends string>({ activeTab, setActiveTab, tabs }: TabHeaderProps<T>) => {
    const theme = useTheme();
    // TODO consider re-arranging this so we have more proper association with tab panels etc.
    return (
        <div className="w-full flex justify-around text-xs font-light" role="tablist">
            {tabs.map((tab) => (
                <button
                    key={tab}
                    onClick={() => setActiveTab(tab)}
                    type="button"
                    role="tab"
                    aria-selected={activeTab === tab}
                    className={clsx(
                        "p-2 text-center border-b-2 transition-colors w-full",
                        // activeTab === tab ? "font-bold" : "font-light"
                        "aria-selected:font-bold",
                        "font-light"
                    )}
                    style={{
                        // we should probably avoid this style object, need to review how theme related styles should be applied
                        borderColor: activeTab === tab ? theme.palette.primary.main : theme.palette.divider,
                        color: activeTab === tab
                            ? theme.palette.primary.main
                            : theme.palette.text.primary,
                    }}
                >
                    {tab}
                </button>
            ))}
        </div>
    );
};

export default TabHeader;