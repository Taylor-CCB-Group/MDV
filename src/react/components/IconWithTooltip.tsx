import { IconButton, type IconButtonProps, Tooltip, type TooltipProps } from "@mui/material";
import type React from "react";
import CustomTooltip from "./CustomTooltip";

export type IconWithTooltipProps = {
    children: React.ReactElement;
    tooltipText: string;
    onClick?: IconButtonProps["onClick"];
    href?: string;
    tooltipProps?: TooltipProps;
    iconButtonProps?: IconButtonProps;
};

const IconWithTooltip = ({
    children,
    tooltipText,
    onClick,
    href,
    tooltipProps,
    iconButtonProps,
}: IconWithTooltipProps) => {
    return (
        <CustomTooltip
            tooltipText={tooltipText}
            {...tooltipProps}
        >
            {href ? (
                <IconButton
                    color="inherit"
                    size="medium"
                    component="a"
                    onClick={onClick}
                    href={href}
                    {...iconButtonProps}
                    style={{
                        color: "inherit",
                        textDecoration: "none",
                        ...iconButtonProps?.style,
                    }}
                >
                    {children}
                </IconButton>
            ) : (
                <IconButton color="inherit" size="medium" onClick={onClick} {...iconButtonProps}>
                    {children}
                </IconButton>
            )}
        </CustomTooltip>
    );
};

export default IconWithTooltip;
