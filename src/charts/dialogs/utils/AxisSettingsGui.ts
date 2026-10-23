import { g } from "@/lib/utils";
import type { AxisConfig, AxisConfig2D } from "@/react/scatter_state";

// premature abstraction?
const axisNames = {
    x: "X Axis",
    y: "Y Axis",
    z: "Z Axis",
    ry: "Right Y Axis",
    tx: "Top X Axis",
};
const sizeLabel = (label: string) => {
    // 🙄 this is a bit dumb
    if (label.toLowerCase().includes("y")) return "Width";
    if (label.toLowerCase().includes("x")) return "Height";
    return "Size";
};

function axisGui(axis: AxisConfig, label: string, callback?: () => void, rotateLabels = false) {
    // rotate_labels is not so useful for scatterplots - numbers tend to look ok
    const rotate = rotateLabels
        ? [
              g({
                  type: "check",
                  label: `Rotate ${label} labels`,
                  current_value: axis.rotate_labels,
                  func: (x) => {
                      axis.rotate_labels = x;
                      callback?.();
                  },
              }),
          ]
        : [];
    return [
        ...rotate,
        g({
            type: "slider",
            max: 20,
            min: 4,
            step: 1,
            label: `${label} text size`,
            current_value: axis.tickfont,
            func: (x) => {
                axis.tickfont = x;
                callback?.();
            },
        }),
        g({
            type: "slider",
            label: `${label} ${sizeLabel(label)}`,
            max: 200,
            min: 20,
            current_value: axis.size,
            func: (x) => {
                axis.size = x;
                callback?.();
            },
        }),
    ];
}

export type AxisGuiOptions = {
    /** Show a "Rotate X axis labels" control, for charts with long category labels. */
    rotateXLabels?: boolean;
};

export default function getAxisGuiSpec(
    config: { axis?: AxisConfig2D },
    callback?: () => void,
    { rotateXLabels = false }: AxisGuiOptions = {},
) {
    if (!config.axis) {
        throw new Error("No axis config found");
    }
    const { axis } = config;
    // mobx doesn't like this Object.entries?
    // const arr = Object.entries(axis).flatMap(([key, value]) => axisGui(value, key, callback));
    const arr = [...axisGui(axis.x, axisNames.x, callback, rotateXLabels), ...axisGui(axis.y, axisNames.y, callback)];
    return g({
        type: "folder",
        label: "Axis Controls",
        current_value: arr,
    });
}
