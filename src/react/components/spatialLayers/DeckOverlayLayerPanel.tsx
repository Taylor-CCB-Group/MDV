import { Typography } from "@mui/material";
import type { DeckOverlayId } from "@/react/spatial_layer_stack";
import { DECK_OVERLAY_LABELS } from "@/react/spatial_layer_stack";

export default function DeckOverlayLayerPanel({ deckId }: { deckId: DeckOverlayId }) {
    return (
        <Typography variant="body2" color="text.secondary">
            {DECK_OVERLAY_LABELS[deckId]} visibility is controlled from the accordion header.
            Opacity is not supported for this overlay type.
        </Typography>
    );
}
