import { IconButton, Tooltip } from "@mui/material";
import type React from "react";

export type IconWithTooltipProps = {
    children: React.ReactNode;
    tooltipText: string;
    onClick: () => void;
    tooltipProps?: object;
    iconButtonProps?: object;
};

const IconWithTooltip = ({ children, tooltipText, onClick, tooltipProps, iconButtonProps }: IconWithTooltipProps) => {
    return (
        <Tooltip
            title={tooltipText}
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
            {...tooltipProps}
        >
            <IconButton color="inherit" size="medium" onClick={onClick} {...iconButtonProps}>
                {children}
            </IconButton>
        </Tooltip>
    );
};

export default IconWithTooltip;
