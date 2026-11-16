import { Checkbox, FormControlLabel, MenuItem, Slider, TextField, Typography } from "@mui/material";
import type { RenderStackGroupBlendMode } from "@spatialdata/layers";

import type { CompositingProps } from "@/react/spatialdata/render_stack_control";

const BLEND_MODES: { value: RenderStackGroupBlendMode; label: string }[] = [
    { value: "pass-through", label: "Pass through" },
    { value: "normal", label: "Normal" },
    { value: "multiply", label: "Multiply" },
    { value: "screen", label: "Screen" },
    { value: "additive", label: "Add" },
    { value: "max", label: "Lighten (max)" },
    { value: "min", label: "Darken (min)" },
];

function isBlendMode(value: string): value is RenderStackGroupBlendMode {
    return BLEND_MODES.some((mode) => mode.value === value);
}

/**
 * Opacity, blend mode and isolation for layers drawn together: a group, or the data
 * points. Any blend mode but pass-through draws the layers as one image first, so it
 * implies isolation.
 */
export default function CompositingControls({
    value,
    onChange,
}: {
    value: CompositingProps;
    onChange: (patch: Partial<CompositingProps>) => void;
}) {
    const blended = value.blendMode !== "pass-through";
    return (
        <div className="grid gap-2">
            <div>
                <Typography variant="body2" id="compositing-opacity">
                    Opacity
                </Typography>
                <Slider
                    size="small"
                    min={0}
                    max={1}
                    step={0.01}
                    value={value.opacity}
                    aria-labelledby="compositing-opacity"
                    valueLabelDisplay="auto"
                    onChange={(_, next) => {
                        if (typeof next === "number") onChange({ opacity: next });
                    }}
                />
            </div>
            <TextField
                select
                size="small"
                label="Blend mode"
                value={value.blendMode}
                onChange={(e) => {
                    if (isBlendMode(e.target.value)) onChange({ blendMode: e.target.value });
                }}
            >
                {BLEND_MODES.map((mode) => (
                    <MenuItem key={mode.value} value={mode.value}>
                        {mode.label}
                    </MenuItem>
                ))}
            </TextField>
            <FormControlLabel
                control={
                    <Checkbox
                        size="small"
                        checked={blended || value.isolate}
                        disabled={blended}
                        onChange={(e) => onChange({ isolate: e.target.checked })}
                    />
                }
                label={<Typography variant="body2">Draw as one layer (isolate)</Typography>}
            />
            <Typography variant="caption" color="text.secondary">
                {blended
                    ? "A blend mode draws these layers as one image, then blends it with what is below."
                    : value.isolate
                      ? "Drawn as one image first, so where layers overlap they don't add up; the opacity applies once."
                      : "Each layer draws on its own, with this opacity applied to it; overlaps add up."}
            </Typography>
        </div>
    );
}
