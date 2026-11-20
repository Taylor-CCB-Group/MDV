import { Checkbox, FormControlLabel, ListSubheader, MenuItem, Slider, TextField, Typography } from "@mui/material";
import type { RenderStackBlendMode, RenderStackGroupBlendMode } from "@spatialdata/layers";
import { useId } from "react";

import {
    BLEND_MODE_LABELS,
    blendModeEffect,
    COLOR_BLEND_MODES,
    COMPOSITING_OPERATORS,
    isCompositingOperator,
    isGroupBlendMode,
    isLeafBlendMode,
} from "@/react/spatialdata/blend_modes";
import type { CompositingProps } from "@/react/spatialdata/render_stack_control";

const NOT_ISOLATED =
    "Drawn as Normal here: it acts on what is beneath it inside a group that is drawn as one layer (one with a blend mode, or isolated), and this is not in one.";

/**
 * The options of a blend mode select: pass-through for groups, the colour modes, and the
 * compositing operators under a heading of their own. A MUI select takes its options as
 * direct children, so this returns an array rather than a fragment.
 */
function blendModeOptions(modes: readonly RenderStackGroupBlendMode[], operators: boolean) {
    const items = modes.map((mode) => (
        <MenuItem key={mode} value={mode}>
            {BLEND_MODE_LABELS[mode]}
        </MenuItem>
    ));
    if (!operators) return items;
    return [
        ...items,
        <ListSubheader key="operators">Masks (inside an isolated group)</ListSubheader>,
        ...COMPOSITING_OPERATORS.map((mode) => (
            <MenuItem key={mode} value={mode}>
                {BLEND_MODE_LABELS[mode]}
            </MenuItem>
        )),
    ];
}

/**
 * Opacity, blend mode and isolation for layers drawn together: a group, or the data
 * points. Any blend mode but pass-through draws the layers as one image first, so it
 * implies isolation.
 */
export default function CompositingControls({
    value,
    onChange,
    insideIsolatedGroup,
}: {
    value: CompositingProps;
    onChange: (patch: Partial<CompositingProps>) => void;
    /**
     * Whether these layers sit inside an isolated group, where the compositing operators act.
     * Leave out where they can never be (the data points), and the operators are not offered.
     */
    insideIsolatedGroup?: boolean;
}) {
    const blended = value.blendMode !== "pass-through";
    const operator = isCompositingOperator(value.blendMode);
    const offerOperators = insideIsolatedGroup !== undefined || operator;
    // several layer dialogs can be open in one document
    const opacityLabelId = useId();
    return (
        <div className="grid gap-2">
            <div>
                <Typography variant="body2" id={opacityLabelId}>
                    Opacity
                </Typography>
                <Slider
                    size="small"
                    min={0}
                    max={1}
                    step={0.01}
                    value={value.opacity}
                    aria-labelledby={opacityLabelId}
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
                    if (isGroupBlendMode(e.target.value)) onChange({ blendMode: e.target.value });
                }}
            >
                {blendModeOptions(["pass-through", ...COLOR_BLEND_MODES], offerOperators)}
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
            <Typography variant="caption" color={operator && !insideIsolatedGroup ? "warning" : "text.secondary"}>
                {operator
                    ? insideIsolatedGroup
                        ? `Drawn as one image, then composited in the group around it: ${blendModeEffect(value.blendMode)}.`
                        : NOT_ISOLATED
                    : blended
                      ? "A blend mode draws these layers as one image, then blends it with what is below."
                      : value.isolate
                        ? "Drawn as one image first, so where layers overlap they don't add up; the opacity applies once."
                        : "Each layer draws on its own, with this opacity applied to it; overlaps add up."}
            </Typography>
        </div>
    );
}

/**
 * A layer's own blend mode: how it blends with the layers beneath it in its group. It acts
 * only inside an isolated group; elsewhere the layer draws normally, and this says so.
 */
export function LayerBlendModeControl({
    value,
    onChange,
    insideIsolatedGroup,
}: {
    value: RenderStackBlendMode;
    onChange: (mode: RenderStackBlendMode) => void;
    insideIsolatedGroup: boolean;
}) {
    const inactive = !insideIsolatedGroup && value !== "normal";
    return (
        <div className="grid gap-1">
            <TextField
                select
                size="small"
                label="Blend mode"
                value={value}
                onChange={(e) => {
                    if (isLeafBlendMode(e.target.value)) onChange(e.target.value);
                }}
            >
                {blendModeOptions(COLOR_BLEND_MODES, true)}
            </TextField>
            <Typography variant="caption" color={inactive ? "warning" : "text.secondary"}>
                {inactive
                    ? NOT_ISOLATED
                    : !insideIsolatedGroup
                      ? "Takes effect inside a group drawn as one layer: one with a blend mode, or isolated. Put this layer in such a group to blend it with the layers beneath it there."
                      : value === "normal"
                        ? "Drawn over the layers beneath it in its group."
                        : `In its group: ${blendModeEffect(value)}.`}
            </Typography>
        </div>
    );
}
