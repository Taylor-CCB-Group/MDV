import { Tooltip, type TooltipProps } from "@mui/material";


export type CustomTooltipProps = {
    tooltipText: string;
    children: React.ReactElement;
    tooltipProps?: TooltipProps; 
};

const CustomTooltip = ({ tooltipText, tooltipProps, children }: CustomTooltipProps) => {
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
            {children}
        </Tooltip>
    );
};

export default CustomTooltip;