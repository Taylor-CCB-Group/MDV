import { Typography } from "@mui/material";
import { findRenderStackEntry } from "@spatialdata/layers";
import { runInAction } from "mobx";
import { observer } from "mobx-react-lite";

import { isInsideIsolatedGroup, leafBlendMode, setLeafBlendMode } from "@/react/spatialdata/blend_modes";
import { useConfig } from "../../hooks";
import type { SpatialDataMdvReactConfig } from "../SpatialDataMDVReact";
import { LayerBlendModeControl } from "./CompositingControls";

/**
 * A SpatialData layer's own blend mode. Written in place with no generation bump: only the
 * stack's structure reads it, and the canvas re-renders on that (`useRenderStackStructure`),
 * so the layer's config and data are untouched.
 */
const LayerBlendModePanel = observer(function LayerBlendModePanel({ entryId }: { entryId: string }) {
    const config = useConfig<SpatialDataMdvReactConfig>();
    const stack = config.renderStack;
    const entry = stack && findRenderStackEntry(stack, entryId);
    if (!stack || entry?.kind !== "spatial") return null;
    return (
        <div className="grid gap-1 mt-3">
            <Typography variant="subtitle2">Compositing</Typography>
            <LayerBlendModeControl
                value={leafBlendMode(entry)}
                insideIsolatedGroup={isInsideIsolatedGroup(stack, entryId)}
                onChange={(mode) => {
                    runInAction(() => {
                        if (config.renderStack) setLeafBlendMode(config.renderStack, entryId, mode);
                    });
                }}
            />
        </div>
    );
});

export default LayerBlendModePanel;
