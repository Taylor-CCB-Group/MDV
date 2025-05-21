import { useTheme } from "@mui/material";
import clsx from "clsx";

type TabHeaderProps<T extends string> = {
    activeTab: T;
    setActiveTab: (tab: T) => void;
    tabs: T[];
};

const TabHeader = <T extends string>({ activeTab, setActiveTab, tabs }: TabHeaderProps<T>) => {
    const theme = useTheme();

    return (
        <div className="w-full flex justify-around text-xs font-light">
            {tabs.map((tab) => (
                <button
                    key={tab}
                    onClick={() => setActiveTab(tab)}
                    type="button"
                    className={clsx(
                        "p-2 text-center border-b-2 transition-colors w-full",
                        activeTab === tab ? "font-bold" : "font-light"
                    )}
                    style={{
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