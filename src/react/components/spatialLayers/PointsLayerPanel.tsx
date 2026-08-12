import { Divider, MenuItem, Slider, TextField, Typography } from "@mui/material";
import { DEFAULT_POINTS_MEMORY_CAP } from "@spatialdata/core";
import { PointsFeatureStateProvider, usePointsFeatureState } from "@spatialdata/vis";
import { observer } from "mobx-react-lite";

import type { PointsLayerConfig, PointsLayerUpdate } from "@/react/spatialdata/points_layer_config";
import { useChart } from "../../context";
import type { SpatialDataMdvReact, SpatialDataMdvReactConfig } from "../SpatialDataMDVReact";
import PointsFeatureFilterPanel from "./PointsFeatureFilterPanel";

type Props = {
    config: PointsLayerConfig;
    updateLayer: PointsLayerUpdate;
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
    return Array.from(new Set([1, 2, 4, 8, 16].map((m) => m * 1_000_000).concat(currentCap))).sort((a, b) => a - b);
}

/**
 * How much of the element is actually in memory, straight from the engine — the
 * context for the memory cap above it, and the only place the user finds out that
 * what they are looking at is capped.
 */
function PointsInMemory({ config }: { config: PointsLayerConfig }) {
    // Engine-backed read that updates on notify; see PointsFeatureFilterPanel.
    "use no memo";
    const { truncation } = usePointsFeatureState(config);
    if (!truncation) return null;
    // Deliberately the size of the batch held in memory, not a per-selection matched
    // count: `loaded` is the covered batch, which overstates a selection that filters
    // that batch in memory. A precise per-selection count needs engine support.
    const message = truncation.truncated
        ? `${truncation.loaded.toLocaleString()}${
              truncation.total !== undefined ? ` of ${truncation.total.toLocaleString()}` : ""
          } points in memory — capped; raise the cap for more.`
        : truncation.filtered
          ? `${truncation.loaded.toLocaleString()} points in memory; view filtered to selection.`
          : `All ${truncation.loaded.toLocaleString()} points loaded (not capped).`;
    return (
        <Typography variant="caption" color={truncation.truncated ? "warning.main" : "text.secondary"}>
            {message}
        </Typography>
    );
}

/**
 * Says out loud that the flat colour above is not what the user is looking at. An
 * element with a feature catalog is coloured per feature (on by default upstream and
 * not switchable — see the note at the colour control), so the swatch only shows
 * through for points with no feature code, while its alpha applies either way.
 */
function FlatColourNote({ config }: { config: PointsLayerConfig }) {
    // Engine-backed read that updates on notify; see PointsFeatureFilterPanel.
    "use no memo";
    const { catalog } = usePointsFeatureState(config);
    if (!catalog) return null;
    return (
        <Typography variant="caption" color="text.secondary">
            Points are coloured per feature, so this sets opacity and the fallback colour for points with no feature.
            Per-feature colours are below.
        </Typography>
    );
}

function PointsStyleControls({ config, updateLayer, engineAvailable }: Props & { engineAvailable: boolean }) {
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

            {/* No colour-by-feature switch. It looks like the obvious companion to the
                per-feature swatches below, but `colorByFeature: false` does not reach the
                deck layer upstream — `useLayerData` spreads the flag only when truthy, so
                an explicit `false` arrives as `undefined` and the shader's
                `!== false` guard keeps colouring. A switch here would be inert.
                See https://github.com/Taylor-CCB-Group/SpatialData.js/issues/147 */}
            <LabeledControl label="Colour">
                <div className="grid gap-1">
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
                    {engineAvailable ? <FlatColourNote config={config} /> : null}
                </div>
            </LabeledControl>

            <LabeledControl label="Memory cap">
                <div className="grid gap-1">
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
                    {engineAvailable ? <PointsInMemory config={config} /> : null}
                </div>
            </LabeledControl>
        </div>
    );
}

/**
 * The layer dialog is a portal with its own React tree, so it cannot see the
 * renderer hook's result. The chart publishes the live engine on
 * `pointsLayerRegistry` and this bridges it back into a provider — the same
 * arrangement as {@link ImageLayerPanel}, and for the same reason: the engine owns
 * the resident window, the feature catalog and the in-flight scans, so a second
 * instance would reload everything and answer different questions from the one
 * that is drawing.
 *
 * The styling controls work without the engine, so they render either way; only
 * the feature filter and the in-memory readout wait for it.
 */
const PointsLayerPanel = observer(function PointsLayerPanel({ config, updateLayer }: Props) {
    const chart = useChart<SpatialDataMdvReactConfig, SpatialDataMdvReact>();
    const registry = chart.pointsLayerRegistry;

    // `engineAvailable` gates the two engine-backed readouts inside the style controls:
    // `usePointsFeatureState` throws outside the provider, so they must not render here.
    if (!registry) {
        return (
            <div className="grid gap-3">
                <PointsStyleControls config={config} updateLayer={updateLayer} engineAvailable={false} />
                <Typography variant="caption" color="text.secondary">
                    Waiting for points data…
                </Typography>
            </div>
        );
    }

    return (
        <PointsFeatureStateProvider engine={registry.engine} target={registry.resolveTarget(config.id)}>
            <div className="grid gap-3">
                <PointsStyleControls config={config} updateLayer={updateLayer} engineAvailable={true} />
                <Divider />
                <PointsFeatureFilterPanel config={config} updateLayer={updateLayer} />
            </div>
        </PointsFeatureStateProvider>
    );
});

export default PointsLayerPanel;
