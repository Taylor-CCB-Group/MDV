/**
 * Feature-row state classification for the points feature filter panel.
 *
 * TEMPORARY MIRROR of `packages/vis/src/SpatialCanvas/featureRowState.ts` upstream.
 * `usePointsFeatureState` (exported since @spatialdata/vis 0.7.0) hands back the raw
 * engine signals — resident, rendered, scanning, whether the element supports an
 * on-demand scan — but the function that *reads* them is not exported from the
 * package entry in 0.7.0, and the exports map has only `"."`, so there is no
 * deep-import route either.
 *
 * Copied rather than re-derived deliberately: the precedence is not obvious (points
 * in memory beat selection and scan state; a deselected-but-loaded feature is cached,
 * not dropped; an unknown resident set has to read as *shown*), and getting it wrong
 * greys a feature the canvas is currently drawing.
 *
 * DELETE THIS FILE when MDV's pin moves past the release carrying
 * https://github.com/Taylor-CCB-Group/SpatialData.js/pull/146 and import
 * `describeFeatureRowState` / `featureRowOpacity` from `@spatialdata/vis` instead.
 * Until then, any edit here is a divergence from what the canvas actually does.
 */

/** Why a feature row is (or isn't) greyed — drives both the dimming and the
 * diagnostic tooltip so they can never disagree. */
export type FeatureRowTone = "resident" | "loaded" | "cached" | "loading" | "noIndex" | "notLoaded";

export interface FeatureRowState {
    tone: FeatureRowTone;
    /** Whether the row is dimmed (its points are not on screen). */
    greyed: boolean;
    /** Short state label, e.g. "loaded", "loading", "not loaded". */
    label: string;
    /** One sentence explaining the state / why it is greyed. */
    reason: string;
}

export interface FeatureRowStateInput {
    /** In the preloaded (resident) window. */
    resident: boolean;
    /** On screen now via the last-completed feature-index scan. */
    rendered: boolean;
    /** In the current selection (checked). */
    selected: boolean;
    /** A feature-index scan for the current selection is in flight. */
    scanning: boolean;
    /** The element can fetch non-resident features on demand (has a feature index). */
    supportsOnDemandLoad: boolean;
    /** The resident set is known (false → we can't distinguish, treat as shown). */
    residentKnown: boolean;
}

/**
 * Classify a feature's render state from the signals the panel already has.
 * Precedence matters: `resident`/`rendered` (its points are in memory) win over
 * selection/scan state. `rendered` here means "in the loaded matched batch",
 * i.e. in memory — a deselected-but-loaded feature is `cached`, not dropped,
 * because removing a feature filters the in-memory batch rather than re-scanning
 * (re-adding it is instant).
 */
export function describeFeatureRowState({
    resident,
    rendered,
    selected,
    scanning,
    supportsOnDemandLoad,
    residentKnown,
}: FeatureRowStateInput): FeatureRowState {
    if (!residentKnown) {
        return {
            tone: "loaded",
            greyed: false,
            label: "shown",
            reason: "The resident set is unknown for this element, so every feature is treated as shown.",
        };
    }
    if (resident) {
        return {
            tone: "resident",
            greyed: false,
            label: "resident",
            reason: "In the preloaded window — shown by filtering the in-memory batch (no dataset scan; a large batch can still take a moment to re-filter).",
        };
    }
    if (rendered) {
        return selected
            ? {
                  tone: "loaded",
                  greyed: false,
                  label: "loaded",
                  reason: "On screen via the feature-index scan for the current selection.",
              }
            : {
                  tone: "cached",
                  greyed: false,
                  label: "in memory",
                  reason: "Loaded in the matched batch but hidden (deselected); re-adding it is instant, no scan.",
              };
    }
    if (selected && scanning) {
        return {
            tone: "loading",
            greyed: true,
            label: "loading",
            reason: "Selected — its feature-index scan is in progress.",
        };
    }
    if (!supportsOnDemandLoad) {
        return {
            tone: "noIndex",
            greyed: true,
            label: "not in sample",
            reason: "Beyond the resident window, and this dataset has no feature index, so it can't be fetched on demand. Raise the memory cap or rewrite the dataset with an index.",
        };
    }
    return {
        tone: "notLoaded",
        greyed: true,
        label: "not loaded",
        reason: "Beyond the resident window; select it to fetch its points via the feature-index scan.",
    };
}

/** Opacity for a row given its state: crisp when its points are on screen,
 * mid-dim while loading, fully dim when not loaded. */
export function featureRowOpacity(state: FeatureRowState): number {
    if (!state.greyed) {
        return 1;
    }
    return state.tone === "loading" ? 0.6 : 0.4;
}
