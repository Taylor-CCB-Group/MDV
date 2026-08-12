import { MenuItem, Slider, TextField, Typography } from "@mui/material";
import { DEFAULT_POINTS_MEMORY_CAP } from "@spatialdata/core";
import type { LayerConfig } from "@spatialdata/vis";

type PointsLayerConfig = Extract<LayerConfig, { type: "points" }>;

type Props = {
    config: PointsLayerConfig;
    updateLayer: (updates: Partial<PointsLayerConfig>) => void;
};

const toHex = (value: [number, number, number, number]) =>
    `#${value
        .slice(0, 3)
        .map((channel) => Math.max(0, Math.min(255, channel)).toString(16).padStart(2, "0"))
        .join("")}`;

const fromHex = (hex: string): [number, number, number] => [
    Number.parseInt(hex.slice(1, 3), 16),
    Number.parseInt(hex.slice(3, 5), 16),
    Number.parseInt(hex.slice(5, 7), 16),
];

function LabeledControl({ label, children }: { label: string; children: React.ReactNode }) {
    return (
        <div className="grid grid-cols-[7rem_minmax(0,1fr)] items-center gap-2">
            <Typography
                fontSize="small"
                sx={{ alignSelf: "center", justifySelf: "end", textAlign: "right", paddingRight: 2 }}
            >
                {label}
            </Typography>
            <div className="min-w-0">{children}</div>
        </div>
    );
}

/**
 * Discrete options rather than a free number input: each choice reloads the
 * resident window, and a text field would do that on every keystroke. The current
 * value is always included so a saved config off the preset list still shows.
 */
function memoryCapOptions(currentCap: number): number[] {
    return Array.from(new Set([1, 2, 4, 8, 16].map((m) => m * 1_000_000).concat(currentCap))).sort(
        (a, b) => a - b,
    );
}

export default function PointsLayerPanel({ config, updateLayer }: Props) {
    const pointSize = config.pointSize ?? 4;
    const color = config.color ?? [100, 149, 237, 200];
    const currentCap = config.pointsMemoryCap ?? DEFAULT_POINTS_MEMORY_CAP;

    return (
        <div className="grid gap-3">
            <LabeledControl label="Point size">
                <div className="flex items-center gap-3">
                    <Slider
                        size="small"
                        min={0.01}
                        max={12}
                        step={0.01}
                        value={pointSize}
                        onChange={(_, value) => {
                            if (typeof value === "number") updateLayer({ pointSize: value });
                        }}
                    />
                    <Typography variant="caption" sx={{ minWidth: "2.5rem" }}>
                        {pointSize.toFixed(2)}
                    </Typography>
                </div>
            </LabeledControl>

            <LabeledControl label="Colour">
                <div className="flex items-center gap-3">
                    <input
                        type="color"
                        value={toHex(color)}
                        onChange={(event) => {
                            const [r, g, b] = fromHex(event.target.value);
                            updateLayer({ color: [r, g, b, color[3] ?? 255] });
                        }}
                    />
                    <Slider
                        size="small"
                        min={0}
                        max={255}
                        step={1}
                        value={color[3] ?? 255}
                        onChange={(_, value) => {
                            if (typeof value === "number") {
                                updateLayer({ color: [color[0], color[1], color[2], value] });
                            }
                        }}
                    />
                </div>
            </LabeledControl>

            <LabeledControl label="Memory cap">
                <TextField
                    select
                    size="small"
                    fullWidth
                    value={currentCap}
                    onChange={(event) => updateLayer({ pointsMemoryCap: Number(event.target.value) })}
                    helperText="Max rows kept in memory. Higher draws more points; picking is limited to ~16.7M per layer."
                >
                    {memoryCapOptions(currentCap).map((cap) => (
                        <MenuItem key={cap} value={cap}>
                            {`${(cap / 1_000_000).toLocaleString(undefined, { maximumFractionDigits: 1 })}M rows`}
                            {cap === DEFAULT_POINTS_MEMORY_CAP ? " (default)" : ""}
                        </MenuItem>
                    ))}
                </TextField>
            </LabeledControl>
        </div>
    );
}
