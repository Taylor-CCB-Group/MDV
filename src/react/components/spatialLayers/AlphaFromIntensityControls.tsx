import { Checkbox, FormControlLabel, Slider, Typography } from "@mui/material";
import { type AlphaFromIntensity, resolveAlphaFromIntensity } from "@spatialdata/layers";
import { useId } from "react";

/** sd.js's own defaults for the option (its `ImageChannelPanel`). */
const DEFAULT_ALPHA_FROM_INTENSITY = { threshold: 0, gain: 1 } as const satisfies AlphaFromIntensity;

/** Gain is set on a log scale, 0.1 to 100, as in sd.js's panel. */
const GAIN_LOG_RANGE = { min: -1, max: 2 };

/**
 * "Opacity from brightness" for an image layer: dark pixels become transparent, so a
 * dark-background image (fluorescence, IMC) can be drawn over another (brightfield, H&E).
 * Writes straight to the entry (`alphaFromIntensity`) on every change, as the other image
 * controls do: the canvas reads it from there.
 */
export default function AlphaFromIntensityControls({
    value,
    onChange,
}: {
    value: unknown;
    onChange: (next: AlphaFromIntensity | undefined) => void;
}) {
    const resolved = resolveAlphaFromIntensity(value);
    const thresholdId = useId();
    const gainId = useId();
    return (
        <div className="grid gap-1">
            <FormControlLabel
                control={
                    <Checkbox
                        size="small"
                        checked={resolved !== null}
                        onChange={(e) => onChange(e.target.checked ? { ...DEFAULT_ALPHA_FROM_INTENSITY } : undefined)}
                    />
                }
                label={<Typography variant="body2">Opacity from brightness</Typography>}
                title="Dark areas become transparent, so this image can be drawn over another, such as H&E"
            />
            {resolved && (
                <>
                    <Typography variant="body2" id={thresholdId}>
                        Threshold
                    </Typography>
                    <Slider
                        size="small"
                        min={0}
                        max={1}
                        step={0.01}
                        value={resolved.threshold}
                        aria-labelledby={thresholdId}
                        valueLabelDisplay="auto"
                        onChange={(_, next) => {
                            if (typeof next === "number") onChange({ ...resolved, threshold: next });
                        }}
                    />
                    <Typography variant="body2" id={gainId}>
                        Gain
                    </Typography>
                    <Slider
                        size="small"
                        min={GAIN_LOG_RANGE.min}
                        max={GAIN_LOG_RANGE.max}
                        step={0.01}
                        value={Math.log10(Math.max(resolved.gain, 10 ** GAIN_LOG_RANGE.min))}
                        scale={(x) => Number((10 ** x).toPrecision(3))}
                        aria-labelledby={gainId}
                        valueLabelDisplay="auto"
                        onChange={(_, next) => {
                            if (typeof next === "number")
                                onChange({ ...resolved, gain: Number((10 ** next).toPrecision(3)) });
                        }}
                    />
                    <Typography variant="caption" color="text.secondary">
                        At or below the threshold the image is transparent; the gain sets how fast it becomes opaque
                        above it. A low threshold with a high gain hides near-black and leaves the rest as it was.
                    </Typography>
                </>
            )}
        </div>
    );
}
