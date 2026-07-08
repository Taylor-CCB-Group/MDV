import { Typography } from "@mui/material";
import type { DeckOverlayId } from "@/react/spatialdata/host_overlay_ids";
import { DECK_OVERLAY_LABELS } from "@/react/spatialdata/host_overlay_ids";

export default function DeckOverlayLayerPanel({ deckId }: { deckId: DeckOverlayId }) {
    return (
        <Typography variant="body2" color="text.secondary">
            {DECK_OVERLAY_LABELS[deckId]} visibility is controlled from the accordion header.
            Opacity is not supported for this overlay type.
        </Typography>
    );
}
