import { Button, Typography } from "@mui/material";
import { observer } from "mobx-react-lite";

import { useRenderStackEntry } from "@/react/spatialdata/render_stack_control";

/**
 * A layer group: what it holds, and taking its layers back out of it. How it is blended is
 * set from its row in the list (`CompositingButton`).
 */
const GroupLayerPanel = observer(function GroupLayerPanel({
    entryId,
    onUngroup,
}: {
    entryId: string;
    onUngroup: (entryId: string) => void;
}) {
    const { entry } = useRenderStackEntry(entryId);
    if (entry?.kind !== "group") return null;
    return (
        <div className="grid gap-2">
            <Typography variant="caption" color="text.secondary">
                {entry.children.length === 0
                    ? "Empty. Drag layers onto it in the list above to add them."
                    : `${entry.children.length} ${entry.children.length === 1 ? "layer" : "layers"}. Drag layers in or out in the list above; drag sideways to change how deep they sit. Set its blending from the icon on its row.`}
            </Typography>
            <div>
                <Button size="small" variant="outlined" onClick={() => onUngroup(entryId)}>
                    Ungroup
                </Button>
            </div>
        </div>
    );
});

export default GroupLayerPanel;
