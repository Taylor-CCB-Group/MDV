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

/**
 * The only config the feature filter reads: the selection, in both its durable and
 * legacy forms, plus the per-feature colour overrides.
 *
 * It is a narrowed type rather than the whole layer config for a performance reason
 * that a comment could not enforce. The layer dialog is an `observer` and
 * `touchRenderStackEntry` subscribes it to *every* key of `entry.props`, so an opacity
 * drag re-renders it; `spatialEntryAsLayerConfig` then hands down a freshly spread
 * object, so nothing downstream can tell the cosmetic edit from a real one. Given the
 * full config the filter panel re-rendered all 541 feature rows per drag step —
 * ~350-450ms of MUI, against ~10ms for the deck update the drag was actually for.
 *
 * Narrowing makes that structural: `memo` can compare the three fields that matter,
 * and TypeScript stops a future edit from quietly reading `pointSize` and putting the
 * dependency back.
 */
export type PointsFeatureFilterConfig = Pick<
    PointsLayerConfig,
    "featureNames" | "featureCodes" | "featureColorOverrides"
>;

/** Same-value check for {@link PointsFeatureFilterConfig}, used as `memo`'s comparator. */
export function samePointsFeatureFilterConfig(a: PointsFeatureFilterConfig, b: PointsFeatureFilterConfig): boolean {
    return (
        sameStringList(a.featureNames, b.featureNames) &&
        sameNumberList(a.featureCodes, b.featureCodes) &&
        sameOverrides(a.featureColorOverrides, b.featureColorOverrides)
    );
}

function sameStringList(a: readonly string[] | undefined, b: readonly string[] | undefined) {
    if (a === b) return true;
    if (!a || !b || a.length !== b.length) return false;
    return a.every((value, index) => value === b[index]);
}

function sameNumberList(a: readonly number[] | undefined, b: readonly number[] | undefined) {
    if (a === b) return true;
    if (!a || !b || a.length !== b.length) return false;
    return a.every((value, index) => value === b[index]);
}

function sameOverrides(
    a: Record<string, [number, number, number]> | undefined,
    b: Record<string, [number, number, number]> | undefined,
) {
    if (a === b) return true;
    if (!a || !b) return false;
    const keys = Object.keys(a);
    if (keys.length !== Object.keys(b).length) return false;
    return keys.every((key) => {
        const left = a[key];
        const right = b[key];
        return right !== undefined && left[0] === right[0] && left[1] === right[1] && left[2] === right[2];
    });
}
