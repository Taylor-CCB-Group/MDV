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

/** Why a layer's mode or a mask does nothing where it is: said the same way everywhere. */
export const ONLY_IN_ONE_LAYER_GROUP = "Only works inside a group drawn as one layer";

/**
 * The options of a blend mode select: pass-through for groups, the colour modes, and the
 * masks under a heading of their own. A MUI select takes its options as direct children, so
 * this returns an array rather than a fragment.
 */
function blendModeOptions(modes: readonly RenderStackGroupBlendMode[], masks: boolean) {
    const items = modes.map((mode) => (
        <MenuItem key={mode} value={mode} dense>
            {BLEND_MODE_LABELS[mode]}
        </MenuItem>
    ));
    if (!masks) return items;
    return [
        ...items,
        <ListSubheader key="masks">Masks</ListSubheader>,
        ...COMPOSITING_OPERATORS.map((mode) => (
            <MenuItem key={mode} value={mode} dense>
                {BLEND_MODE_LABELS[mode]}
            </MenuItem>
        )),
    ];
}

function Hint({ children, warn = false }: { children: string; warn?: boolean }) {
    return (
        <Typography variant="caption" color={warn ? "warning.main" : "text.secondary"} sx={{ display: "block" }}>
            {children}
        </Typography>
    );
}

/** One plain line on what a group's settings do, or why they do nothing where it is. */
export function groupCompositingHint(value: CompositingProps, insideOneLayerGroup: boolean | undefined) {
    const mode = value.blendMode;
    if (isCompositingOperator(mode)) {
        return insideOneLayerGroup
            ? { text: `Drawn as one layer, then masks what is below it in the outer group: ${blendModeEffect(mode)}.` }
            : { text: `${ONLY_IN_ONE_LAYER_GROUP}; here it draws as Normal.`, warn: true };
    }
    if (mode !== "pass-through") return { text: `Drawn as one layer, then ${blendModeEffect(mode)}.` };
    return value.isolate
        ? { text: "Drawn as one layer, then over what is below it, so overlaps inside do not add up." }
        : { text: "Each layer blends with what is below it on its own." };
}

/**
 * How a group, or the data points, are drawn over what is below them: blend mode, and whether
 * they are drawn as one layer first (implied by any mode but pass-through). Opacity only where
 * the row has no slider for it (the data points).
 */
export default function CompositingControls({
    value,
    onChange,
    insideOneLayerGroup,
    showOpacity = false,
}: {
    value: CompositingProps;
    onChange: (patch: Partial<CompositingProps>) => void;
    /** Whether masks act here. Leave out where they never can (the data points) to not offer them. */
    insideOneLayerGroup?: boolean;
    showOpacity?: boolean;
}) {
    const blended = value.blendMode !== "pass-through";
    const masks = insideOneLayerGroup !== undefined || isCompositingOperator(value.blendMode);
    const opacityLabelId = useId();
    const hint = groupCompositingHint(value, insideOneLayerGroup);
    return (
        <div className="grid gap-1">
            {showOpacity && (
                <div>
                    <Typography variant="caption" id={opacityLabelId}>
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
            )}
            <TextField
                select
                size="small"
                label="Blend mode"
                value={value.blendMode}
                onChange={(e) => {
                    if (isGroupBlendMode(e.target.value)) onChange({ blendMode: e.target.value });
                }}
            >
                {blendModeOptions(["pass-through", ...COLOR_BLEND_MODES], masks)}
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
                label={<Typography variant="body2">Draw as one layer</Typography>}
            />
            <Hint warn={hint.warn}>{hint.text}</Hint>
        </div>
    );
}

/** One plain line on what a layer's own mode does, or why it does nothing where it is. */
export function layerBlendHint(value: RenderStackBlendMode, insideOneLayerGroup: boolean) {
    if (!insideOneLayerGroup) {
        return value === "normal"
            ? { text: `${ONLY_IN_ONE_LAYER_GROUP}.` }
            : { text: `${ONLY_IN_ONE_LAYER_GROUP}; here it draws as Normal.`, warn: true };
    }
    if (value === "normal") return { text: "Drawn over the layers below it in its group." };
    return { text: `In its group, ${blendModeEffect(value)}.` };
}

/** A layer's own blend mode: how it blends with the layers below it in its group. */
export function LayerBlendModeControl({
    value,
    onChange,
    insideOneLayerGroup,
}: {
    value: RenderStackBlendMode;
    onChange: (mode: RenderStackBlendMode) => void;
    insideOneLayerGroup: boolean;
}) {
    const hint = layerBlendHint(value, insideOneLayerGroup);
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
            <Hint warn={hint.warn}>{hint.text}</Hint>
        </div>
    );
}
