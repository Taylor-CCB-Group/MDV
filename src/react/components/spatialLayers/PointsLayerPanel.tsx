import { Slider, TextField } from "@mui/material";
import type { LayerConfig } from "@spatialdata/vis";

type PointsLayerConfig = Extract<LayerConfig, { type: "points" }>;

type Props = {
    config: PointsLayerConfig;
    updateLayer: (updates: Partial<PointsLayerConfig>) => void;
};

export default function PointsLayerPanel({ config, updateLayer }: Props) {
    const pointSize = config.pointSize ?? 4;
    const color = config.color ?? [100, 149, 237, 200];

    return (
        <div className="grid gap-3">
            <div className="flex items-center gap-3">
                <span className="w-20 text-xs uppercase tracking-wide text-[hsl(var(--muted-foreground))]">
                    Size
                </span>
                <Slider
                    size="small"
                    min={1}
                    max={20}
                    step={0.5}
                    value={pointSize}
                    onChange={(_, value) => {
                        if (typeof value === "number") {
                            updateLayer({ pointSize: value });
                        }
                    }}
                />
            </div>
            <div className="grid grid-cols-4 gap-2">
                {(["R", "G", "B", "A"] as const).map((label, index) => (
                    <TextField
                        key={label}
                        size="small"
                        type="number"
                        label={label}
                        value={color[index] ?? 0}
                        onChange={(event) => {
                            const next = [...color] as [number, number, number, number];
                            next[index] = Number(event.target.value);
                            updateLayer({ color: next });
                        }}
                    />
                ))}
            </div>
        </div>
    );
}
