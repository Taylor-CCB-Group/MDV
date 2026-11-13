import { Typography } from "@mui/material";
import { useMemo } from "react";
import type { DeckOverlayId } from "@/react/spatialdata/host_overlay_ids";
import { useChart } from "../../context";
import { SettingsGuiList } from "../SettingsDialogComponent";
import type { SpatialDataMdvReact, SpatialDataMdvReactConfig } from "../SpatialDataMDVReact";

/**
 * The MDV data points: colour, tooltip, point size, filtering and density. These used to
 * be chart settings; they are edited here so every layer's appearance is in one dialog.
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
        </div>
    );
}

export default function DeckOverlayLayerPanel({ deckId }: { deckId: DeckOverlayId }) {
    if (deckId === "scatter") return <DataPointsLayerPanel />;
    return (
        <Typography variant="body2" color="text.secondary">
            {deckId === "gate_labels"
                ? "Names of the gates drawn on this chart. Gates are created with the selection tools."
                : "Gates drawn on this chart. Gates are created with the selection tools."}
        </Typography>
    );
}
