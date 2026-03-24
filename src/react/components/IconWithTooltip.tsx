import { IconButton, type IconButtonProps, Tooltip, type TooltipProps } from "@mui/material";
import type React from "react";
import CustomTooltip from "./CustomTooltip";

export type IconWithTooltipProps = {
    children: React.ReactElement;
    tooltipText: string;
    onClick: React.MouseEventHandler<HTMLButtonElement>;
    tooltipProps?: TooltipProps;
    iconButtonProps?: IconButtonProps;
};

const IconWithTooltip = ({ children, tooltipText, onClick, tooltipProps, iconButtonProps }: IconWithTooltipProps) => {
    return (
        <CustomTooltip
            tooltipText={tooltipText}
            {...tooltipProps}
        >
            <IconButton color="inherit" size="medium" onClick={onClick} {...iconButtonProps}>
                {children}
            </IconButton>
        </CustomTooltip>
    );
};

export default IconWithTooltip;
