import JoinInnerIcon from "@mui/icons-material/JoinInner";
import { ButtonBase, Popover, Typography } from "@mui/material";
import { findRenderStackEntry } from "@spatialdata/layers";
import { runInAction } from "mobx";
import { observer } from "mobx-react-lite";
import { useState } from "react";

import {
    compositingBadge,
    isInsideIsolatedGroup,
    leafBlendMode,
    setLeafBlendMode,
} from "@/react/spatialdata/blend_modes";
import { deckHostLayerId } from "@/react/spatialdata/host_overlay_ids";
import {
    type CompositingProps,
    DEFAULT_COMPOSITING,
    getDataPointsGroup,
    setDataPointsCompositing,
    useRenderStackEntry,
} from "@/react/spatialdata/render_stack_control";
import { useChart } from "../../context";
import { useConfig } from "../../hooks";
import type { SpatialDataMdvReact, SpatialDataMdvReactConfig } from "../SpatialDataMDVReact";
import CompositingControls, { LayerBlendModeControl } from "./CompositingControls";

const DATA_POINTS_ENTRY_ID = deckHostLayerId("scatter");

/** Whether a layer list row has compositing to set: SpatialData layers, groups and the data points. */
export function hasCompositing(entryKind: string, entryId: string): boolean {
    return entryKind === "spatial" || entryKind === "group" || entryId === DATA_POINTS_ENTRY_ID;
}

/**
 * A layer list row's compositing, set from the row itself: the row's badge (its mode, when not
 * the default) or a faint icon, opening a small popover with the controls. Keeps blending out
 * of the properties pane, so it is not at the bottom of every panel.
 */
export const CompositingButton = observer(function CompositingButton({
    entryId,
    name,
}: {
    entryId: string;
    name: string;
}) {
    const config = useConfig<SpatialDataMdvReactConfig>();
    const [anchor, setAnchor] = useState<HTMLElement | null>(null);
    const stack = config.renderStack;
    const entry = stack && findRenderStackEntry(stack, entryId);
    if (!stack || !entry) return null;
    const badge = compositingBadge(stack, entry);
    const open = (e: React.MouseEvent<HTMLElement>) => {
        e.stopPropagation();
        setAnchor(e.currentTarget);
    };
    const label = `blending of ${name}`;
    return (
        <>
            {/* one element whichever it shows, as the popover is anchored to it while a change of
                mode turns the icon into a badge or back */}
            <ButtonBase
                onClick={open}
                aria-label={label}
                aria-haspopup="dialog"
                title={badge?.title ?? "Blending"}
                sx={{
                    flexShrink: 0,
                    borderRadius: 1,
                    ...(badge ? {} : { p: "5px", opacity: 0.45, "&:hover": { opacity: 1 } }),
                }}
            >
                {badge ? (
                    <Typography
                        variant="caption"
                        className={`rounded border px-1 leading-4 select-none ${badge.inactive ? "border-dashed opacity-60" : "border-current"}`}
                        color={badge.inactive ? "text.secondary" : "primary"}
                    >
                        {badge.label}
                    </Typography>
                ) : (
                    <JoinInnerIcon fontSize="small" />
                )}
            </ButtonBase>
            <Popover
                open={anchor !== null}
                anchorEl={anchor}
                onClose={() => setAnchor(null)}
                onClick={(e) => e.stopPropagation()}
                anchorOrigin={{ vertical: "bottom", horizontal: "right" }}
                transformOrigin={{ vertical: "top", horizontal: "right" }}
                slotProps={{ paper: { sx: { p: 1.5, width: 260 }, role: "dialog", "aria-label": label } }}
            >
                {entryId === DATA_POINTS_ENTRY_ID ? (
                    <DataPointsCompositing />
                ) : entry.kind === "group" ? (
                    <GroupCompositing entryId={entryId} />
                ) : (
                    <LayerCompositing entryId={entryId} />
                )}
            </Popover>
        </>
    );
});

const GroupCompositing = observer(function GroupCompositing({ entryId }: { entryId: string }) {
    const config = useConfig<SpatialDataMdvReactConfig>();
    const { entry, patchProps } = useRenderStackEntry(entryId);
    if (entry?.kind !== "group" || !config.renderStack) return null;
    const value: CompositingProps = {
        opacity: entry.props.opacity ?? DEFAULT_COMPOSITING.opacity,
        blendMode: entry.props.blendMode ?? DEFAULT_COMPOSITING.blendMode,
        isolate: entry.props.isolate ?? DEFAULT_COMPOSITING.isolate,
    };
    return (
        <CompositingControls
            value={value}
            onChange={(patch) => patchProps(patch)}
            insideOneLayerGroup={isInsideIsolatedGroup(config.renderStack, entryId)}
        />
    );
});

/**
 * A SpatialData layer's own blend mode. Written in place with no generation bump: only the
 * stack's structure reads it, and the canvas re-renders on that (`useRenderStackStructure`),
 * so the layer's config and data are untouched.
 */
const LayerCompositing = observer(function LayerCompositing({ entryId }: { entryId: string }) {
    const config = useConfig<SpatialDataMdvReactConfig>();
    const stack = config.renderStack;
    const entry = stack && findRenderStackEntry(stack, entryId);
    if (!stack || entry?.kind !== "spatial") return null;
    return (
        <LayerBlendModeControl
            value={leafBlendMode(entry)}
            insideOneLayerGroup={isInsideIsolatedGroup(stack, entryId)}
            onChange={(mode) => {
                runInAction(() => {
                    if (config.renderStack) setLeafBlendMode(config.renderStack, entryId, mode);
                });
            }}
        />
    );
});

/**
 * The data points (filtered-out and coloured) drawn over what is below them, as one layer.
 * Backed by a group that exists only while this is not the default.
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
    return <CompositingControls value={value} onChange={onChange} showOpacity />;
});
