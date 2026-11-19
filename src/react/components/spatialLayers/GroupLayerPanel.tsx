import { Button, Typography } from "@mui/material";
import { observer } from "mobx-react-lite";

import {
    type CompositingProps,
    DEFAULT_COMPOSITING,
    useRenderStackEntry,
} from "@/react/spatialdata/render_stack_control";
import CompositingControls from "./CompositingControls";

/** A layer group: how its layers are composited, and taking them back out of it. */
const GroupLayerPanel = observer(function GroupLayerPanel({
    entryId,
    onUngroup,
}: {
    entryId: string;
    onUngroup: (entryId: string) => void;
}) {
    const { entry, patchProps } = useRenderStackEntry(entryId);
    if (entry?.kind !== "group") return null;
    const value: CompositingProps = {
        opacity: entry.props.opacity ?? DEFAULT_COMPOSITING.opacity,
        blendMode: entry.props.blendMode ?? DEFAULT_COMPOSITING.blendMode,
        isolate: entry.props.isolate ?? DEFAULT_COMPOSITING.isolate,
    };
    return (
        <div className="grid gap-2">
            <Typography variant="caption" color="text.secondary">
                {entry.children.length === 0
                    ? "Empty. Drag layers onto it in the list above to add them."
                    : `${entry.children.length} ${entry.children.length === 1 ? "layer" : "layers"}. Drag layers in or out in the list above; drag sideways to change how deep they sit.`}
            </Typography>
            <CompositingControls value={value} onChange={(patch) => patchProps(patch)} />
            <div>
                <Button size="small" variant="outlined" onClick={() => onUngroup(entryId)}>
                    Ungroup
                </Button>
            </div>
        </div>
    );
});

export default GroupLayerPanel;
