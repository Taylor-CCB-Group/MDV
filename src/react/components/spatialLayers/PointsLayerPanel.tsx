import ExpandMoreIcon from "@mui/icons-material/ExpandMore";
import {
    Accordion,
    AccordionDetails,
    AccordionSummary,
    Divider,
    MenuItem,
    Slider,
    TextField,
    Typography,
} from "@mui/material";
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

/** Cornflower blue, the flat colour a points layer starts with. */
const DEFAULT_POINT_COLOR: [number, number, number, number] = [100, 149, 237, 200];

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
 * How much of the element is actually in memory, straight from the engine — and the
 * only place the user finds out that what they are looking at is capped.
 *
 * It sits in the Advanced *summary* rather than inside it, so collapsing the section
 * hides the cap control without hiding the fact that a cap is biting.
 */
function PointsInMemory({ config }: { config: PointsLayerConfig }) {
    // Engine-backed read that updates on notify; see PointsFeatureFilterPanel.
    "use no memo";
    const { truncation } = usePointsFeatureState(config);
    if (!truncation) return null;
    // Deliberately the size of the batch held in memory, not a per-selection matched
    // count: `loaded` is the covered batch, which overstates a selection that filters
    // that batch in memory. A precise per-selection count needs engine support.
    const loaded = truncation.loaded.toLocaleString();
    const message = truncation.truncated
        ? truncation.total !== undefined
            ? `${loaded} of ${truncation.total.toLocaleString()} points in memory (capped)`
            : `${loaded} points in memory (capped)`
        : truncation.filtered
          ? `${loaded} points in memory · filtered to selection`
          : `All ${loaded} points loaded`;
    return (
        <Typography variant="caption" color={truncation.truncated ? "warning.main" : "text.secondary"} noWrap>
            {message}
        </Typography>
    );
}

function FlatColourControl({
    label,
    color,
    updateLayer,
    note,
}: {
    label: string;
    color: [number, number, number, number];
    updateLayer: PointsLayerUpdate;
    note?: string;
}) {
    return (
        <LabeledControl label={label}>
            <div className="grid justify-items-start gap-1">
                <input
                    type="color"
                    aria-label={label}
                    value={toHex(color)}
                    onChange={(event) => {
                        const [r, g, b] = fromHex(event.target.value);
                        // Alpha is carried through untouched. Layer-wide opacity is the
                        // render stack's own `opacity` prop, driven by the slider on the
                        // layer row — this panel used to offer a second slider for
                        // `color[3]`, which multiplied against it for no gain.
                        updateLayer({ color: [r, g, b, color[3] ?? 255] });
                    }}
                />
                {note ? (
                    <Typography variant="caption" color="text.secondary">
                        {note}
                    </Typography>
                ) : null}
            </div>
        </LabeledControl>
    );
}

type Placement = "inline" | "advanced";

/**
 * The flat colour, placed by how much it can actually do.
 *
 * With a feature catalog the points take their colours from the feature palette, and
 * that colouring cannot even be switched off (upstream SpatialData.js#147), so this
 * swatch only shows through for points carrying no feature code — near enough inert,
 * and it belongs under Advanced. Without a catalog it is the layer's only colour and
 * has to stay in view.
 *
 * `catalog === null` is the one state that positively means "this element has no
 * catalog"; while the answer is still unknown the swatch sits under Advanced, so the
 * only move it can ever make is the revealing one.
 */
function FlatColour({ config, updateLayer, placement }: Props & { placement: Placement }) {
    // Engine-backed read that updates on notify; see PointsFeatureFilterPanel.
    "use no memo";
    const { catalog } = usePointsFeatureState(config);
    const belongs: Placement = catalog === null ? "inline" : "advanced";
    if (belongs !== placement) return null;
    return (
        <FlatColourControl
            label={placement === "inline" ? "Colour" : "Fallback colour"}
            color={config.color ?? DEFAULT_POINT_COLOR}
            updateLayer={updateLayer}
            note={
                placement === "advanced"
                    ? "Points are coloured per feature; this shows only for points with no feature."
                    : undefined
            }
        />
    );
}

function PointsStyleControls({ config, updateLayer, engineAvailable }: Props & { engineAvailable: boolean }) {
    const pointSize = config.pointSize ?? 4;

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

            {/* Without the engine we cannot ask whether a feature catalog exists, so the
                swatch stays in view — hiding the only colour control of a layer that turns
                out to have no catalog is the worse mistake. */}
            {engineAvailable ? (
                <FlatColour config={config} updateLayer={updateLayer} placement="inline" />
            ) : (
                <FlatColourControl
                    label="Colour"
                    color={config.color ?? DEFAULT_POINT_COLOR}
                    updateLayer={updateLayer}
                />
            )}
        </div>
    );
}

/**
 * Everything that is about the *loading* of this layer rather than the look of it.
 *
 * The memory cap used to sit third in the style block, full width with two lines of
 * explanation under it — more prominence than a control most users never touch, and
 * ahead of the feature list they came for.
 */
function PointsAdvanced({ config, updateLayer, engineAvailable }: Props & { engineAvailable: boolean }) {
    const currentCap = config.pointsMemoryCap ?? DEFAULT_POINTS_MEMORY_CAP;

    return (
        <Accordion
            disableGutters
            elevation={0}
            defaultExpanded={false}
            // Collapsed means gone: the cap select is a load-triggering control and has no
            // business in the tab order while the section is shut.
            slotProps={{ transition: { unmountOnExit: true } }}
            sx={{
                border: "1px solid hsl(var(--border))",
                borderRadius: "6px",
                backgroundColor: "transparent",
                backgroundImage: "none",
                "&:before": { display: "none" },
            }}
        >
            <AccordionSummary
                expandIcon={<ExpandMoreIcon fontSize="small" />}
                sx={{
                    px: 1,
                    minHeight: "unset",
                    "& .MuiAccordionSummary-content": {
                        my: 0.75,
                        minWidth: 0,
                        alignItems: "baseline",
                        gap: 1,
                    },
                    "& .MuiAccordionSummary-expandIconWrapper": { color: "hsl(var(--muted-foreground))" },
                }}
            >
                <Typography variant="caption">Advanced</Typography>
                {engineAvailable ? <PointsInMemory config={config} /> : null}
            </AccordionSummary>
            <AccordionDetails sx={{ px: 1, pt: 0, pb: 1.5 }}>
                <div className="grid gap-3">
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
                    {engineAvailable ? (
                        <FlatColour config={config} updateLayer={updateLayer} placement="advanced" />
                    ) : null}
                </div>
            </AccordionDetails>
        </Accordion>
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

    // `engineAvailable` gates the engine-backed readouts inside the style controls and
    // the Advanced section: `usePointsFeatureState` throws outside the provider, so they
    // must not render here.
    if (!registry) {
        return (
            <div className="grid gap-3">
                <PointsStyleControls config={config} updateLayer={updateLayer} engineAvailable={false} />
                <Typography variant="caption" color="text.secondary">
                    Waiting for points data…
                </Typography>
                <PointsAdvanced config={config} updateLayer={updateLayer} engineAvailable={false} />
            </div>
        );
    }

    return (
        <PointsFeatureStateProvider engine={registry.engine} target={registry.resolveTarget(config.id)}>
            <div className="grid gap-3">
                <PointsStyleControls config={config} updateLayer={updateLayer} engineAvailable={true} />
                <Divider />
                <PointsFeatureFilterPanel config={config} updateLayer={updateLayer} />
                <PointsAdvanced config={config} updateLayer={updateLayer} engineAvailable={true} />
            </div>
        </PointsFeatureStateProvider>
    );
});

export default PointsLayerPanel;
