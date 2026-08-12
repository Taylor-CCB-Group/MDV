import type { useSpatialCanvasRendererFromLayerInputs } from "@spatialdata/vis";

/**
 * Derived from the renderer hook rather than imported by name: `@spatialdata/vis`
 * exports `PointsDataEngine` / `PointsLoadTarget` from `SpatialCanvas/public.ts`,
 * but its package entry does not re-export them (0.6.0), and the exports map has
 * only `"."` so there is no deep-import route either. The renderer's return type
 * carries both structurally, which is enough to name them here.
 */
type SpatialCanvasRenderer = ReturnType<typeof useSpatialCanvasRendererFromLayerInputs>;
type PointsDataEngine = SpatialCanvasRenderer["pointsEngine"];
type PointsLoadTarget = ReturnType<SpatialCanvasRenderer["resolvePointsTarget"]>;

/**
 * The live points engine, published from the chart so the layer dialog can reach it.
 *
 * The dialog is its own React tree — a portal with its own `SpatialDataProvider` —
 * so it cannot see the renderer hook's result, and the engine is not something to
 * rebuild on the other side: it owns the resident window, the feature catalog and
 * the in-flight scans, and a second instance would load everything again and answer
 * different questions from the one that is actually drawing.
 *
 * Same shape as {@link ImageLayerRegistry} for the same reason, and published the
 * same way (`chart.setPointsLayerRegistry` in an effect, cleared on unmount).
 *
 * `engine` is deliberately not `observable` beyond the ref: it is a mutable object
 * that notifies its own subscribers, and `usePointsFeatureState` subscribes to it
 * directly. MobX only needs to know when the registry itself is swapped.
 */
export type PointsLayerRegistry = {
    engine: PointsDataEngine;
    /**
     * Layer id → the engine's load target. Must come from the renderer hook rather
     * than being reconstructed here, so panel reads hit the same cache keys the
     * render path writes.
     */
    resolveTarget: (layerId: string) => PointsLoadTarget;
};

export function createPointsLayerRegistry(
    engine: PointsDataEngine,
    resolveTarget: PointsLayerRegistry["resolveTarget"],
): PointsLayerRegistry {
    return { engine, resolveTarget };
}
