import type { LayerConfig } from "@spatialdata/vis";

/**
 * A points layer's serialisable config — the `points` arm of the upstream
 * `LayerConfig` union.
 *
 * Its own module so the two panels that share it (the styling controls and the
 * feature filter) do not have to import each other for a type, and so neither
 * component file carries a non-component export.
 */
export type PointsLayerConfig = Extract<LayerConfig, { type: "points" }>;

export type PointsLayerUpdate = (updates: Partial<PointsLayerConfig>) => void;
