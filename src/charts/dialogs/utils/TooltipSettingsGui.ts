import { g } from "@/lib/utils";
import type { ScatterPlotConfig } from "@/react/scatter_state";

// this will be for more general tooltip settings rather than just scatterplot
export default function getTooltipSettings(config: ScatterPlotConfig, callback?: () => void) {
    return g({
        type: "folder",
        label: "Tooltip",
        current_value: [
            g({
                type: "check",
                label: "Show tooltip",
                current_value: config.tooltip.show,
                func: (x) => {
                    config.tooltip.show = x;
                    callback?.();
                },
            }),
            g({
                type: "column",
                label: "Tooltip value",
                //@ts-expect-error - need a way of dealing with optional column...
                current_value: config.tooltip.column,
                func: async (x) => {
                    //@ts-expect-error pending tooltip column being FieldSpec(s)
                    config.tooltip.column = x;
                    callback?.();
                }
            }),
        ]
    });
}
