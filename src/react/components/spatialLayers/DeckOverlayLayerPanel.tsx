import { Typography } from "@mui/material";
import { runInAction } from "mobx";
import { observer } from "mobx-react-lite";
import { useMemo } from "react";
import type { DeckOverlayId } from "@/react/spatialdata/host_overlay_ids";
import {
    type CompositingProps,
    DEFAULT_COMPOSITING,
    getDataPointsGroup,
    setDataPointsCompositing,
} from "@/react/spatialdata/render_stack_control";
import { useChart } from "../../context";
import { useConfig } from "../../hooks";
import { SettingsGuiList } from "../SettingsDialogComponent";
import type { SpatialDataMdvReact, SpatialDataMdvReactConfig } from "../SpatialDataMDVReact";
import CompositingControls from "./CompositingControls";

/**
 * The MDV data points: colour, tooltip, point size and filtering. These used to be chart
 * settings; they are edited here so every layer's appearance is in one dialog.
 */
function DataPointsLayerPanel() {
    const chart = useChart<SpatialDataMdvReactConfig, SpatialDataMdvReact>();
    // built once per mount: the controls write into these spec objects
    const specs = useMemo(() => chart.getDataLayerSettings(), [chart]);
    return (
        <div className="grid gap-1">
            <Typography variant="caption" color="text.secondary">
                Rows of '{chart.dataStore.name}', drawn at their positions in this region. Always drawn above
                SpatialData layers.
            </Typography>
            <SettingsGuiList specs={specs} />
            <DataPointsCompositing />
        </div>
    );
}

/**
 * How the data points (filtered-out and coloured) are drawn over what is below them,
 * as one layer. Backed by a group that exists only while this is not the default.
 */
const DataPointsCompositing = observer(function DataPointsCompositing() {
    const chart = useChart<SpatialDataMdvReactConfig, SpatialDataMdvReact>();
    const config = useConfig<SpatialDataMdvReactConfig>();
    const stack = config.renderStack;
    const group = stack && getDataPointsGroup(stack);
    const value: CompositingProps = {
        opacity: group?.props.opacity ?? DEFAULT_COMPOSITING.opacity,
        blendMode: group?.props.blendMode ?? DEFAULT_COMPOSITING.blendMode,
        isolate: group?.props.isolate ?? DEFAULT_COMPOSITING.isolate,
    };
    const onChange = (patch: Partial<CompositingProps>) => {
        runInAction(() => {
            if (!config.renderStack) return;
            const structural = setDataPointsCompositing(config.renderStack, { ...value, ...patch });
            if (structural) chart.bumpRenderStackGeneration();
            else chart.bumpRenderStackPropsGeneration();
        });
    };
    return (
        <div className="grid gap-1 mt-2">
            <Typography variant="subtitle2">Compositing</Typography>
            <CompositingControls value={value} onChange={onChange} />
        </div>
    );
});

/** Density contours of the data points, by category or weighted by numeric columns. */
function DensityLayerPanel() {
    const chart = useChart<SpatialDataMdvReactConfig, SpatialDataMdvReact>();
    const specs = useMemo(() => chart.getDensityLayerSettings(), [chart]);
    return (
        <div className="grid gap-1">
            <Typography variant="caption" color="text.secondary">
                Density of the rows of '{chart.dataStore.name}' that pass the filters: for chosen categories of a
                column, or weighted by numeric columns. Drawn whether or not the data points are shown.
            </Typography>
            <SettingsGuiList specs={specs} />
        </div>
    );
}

export default function DeckOverlayLayerPanel({ deckId }: { deckId: DeckOverlayId }) {
    if (deckId === "scatter") return <DataPointsLayerPanel />;
    if (deckId === "density") return <DensityLayerPanel />;
    return (
        <Typography variant="body2" color="text.secondary">
            {deckId === "gate_labels"
                ? "Names of the gates drawn on this chart. Gates are created with the selection tools."
                : "Gates drawn on this chart. Gates are created with the selection tools."}
        </Typography>
    );
}
