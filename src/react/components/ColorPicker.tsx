// nb, was going to use `react-colorful`, for now <input type="color" /> is ok
// https://codesandbox.io/s/opmco?file=/src/PopoverPicker.js
// may in future want to go more in-depth with color picker component.

import { useMemo } from "react";

type MyRgbColor = [number, number, number];

export const PopoverPicker = ({
    color,
    onChange,
}: { color?: MyRgbColor; onChange?: (newColor: MyRgbColor) => void }) => {
    const [r, g, b] = color || [255, 255, 255];

    const backgroundColor = color
        ? `#${r.toString(16).padStart(2, "0")}${g.toString(16).padStart(2, "0")}${b.toString(16).padStart(2, "0")}`
        : "#000000";
    const style = useMemo(
        () => ({
            // padding: '10px',
            border: "none",
            cursor: "pointer",
            width: "100%",
            backgroundColor,
        }),
        [backgroundColor],
    );
    return (
        <input
            type="color"
            className="color-picker-no-border"
            style={style}
            value={backgroundColor}
            onChange={(e) => {
                const c = e.target.value;
                //get rgb values from hex string
                const [r, g, b] = [
                    c.slice(1, 3),
                    c.slice(3, 5),
                    c.slice(5, 7),
                ].map((x) => Number.parseInt(x, 16));
                onChange([r, g, b]);
            }}
        />
    );
};
