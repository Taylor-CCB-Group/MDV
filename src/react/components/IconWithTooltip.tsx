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
                        color: "black",
                    },
                },
                tooltip: {
                    sx: {
                        backgroundColor: "black",
                        fontSize: "0.8rem",
                    },
                },
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
