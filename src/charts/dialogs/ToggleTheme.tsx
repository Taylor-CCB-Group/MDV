import { Brightness4 as Brightness4Icon, Brightness7 as Brightness7Icon } from "@mui/icons-material";
import { IconButton, Tooltip } from "@mui/material";
import { observer } from "mobx-react-lite";

const ToggleTheme = observer(() => {
    const cm = window.mdv.chartManager;
    const { theme } = cm;

    // todo theme should be `"system" | "dark" | "light"`, and stored in local storage
    const toggleTheme = () => {
        if (theme === "dark") cm.setTheme("light");
        else cm.setTheme("dark");
    };

    return (
        <Tooltip
            title="Toggle Theme"
            arrow
            slotProps={{
                arrow: {
                    sx: {
                        color: "var(--tooltip_background_color)"
                    }
                },
                tooltip: {
                    sx: {
                        backgroundColor: "var(--tooltip_background_color)",
                        color: "var(--tooltip_text_color)",
                        fontSize: "0.8rem",
                        fontWeight: "normal",
                    }
                }
            }}
        >
            <IconButton sx={{ ml: 1 }} onClick={toggleTheme} color="inherit">
                {theme === "dark" ? <Brightness4Icon /> : <Brightness7Icon />}
            </IconButton>
        </Tooltip>
    );
});

const ToggleThemeWrapper = () => {
    return <ToggleTheme />;
};

export default ToggleThemeWrapper;
