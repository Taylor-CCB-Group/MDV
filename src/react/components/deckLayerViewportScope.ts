import { isDensityGridViewport, matchesDensityGridView } from "./densityGridUtils";

/**
 * Declares how a deck layer participates in multi-viewport charts (single scatter vs density grid).
 *
 * - `chart-shared`: one layer instance, same geometry/data in every viewport (points, gates, JSON,
 *   selection). deck.gl `layerFilter` draws it in each cell; avoids N clones and duplicate GPU buffers.
 * - `per-viewport`: viewport-specific visuals (density/contour per field today; per-cell scatter tint later).
 *
 * Set on layer props via {@link tagDeckLayerViewportScope}. Prefer tagging at layer creation over id heuristics.
 */
export const MDV_DECK_LAYER_VIEWPORT_SCOPE = "mdvDeckLayerViewportScope" as const;

export type DeckLayerViewportScope = "chart-shared" | "per-viewport";

export type DeckLayerScopeInput = {
    id: string;
    props?: Record<string, unknown> | null;
};

export function getDeckLayerViewportScope(layer: DeckLayerScopeInput): DeckLayerViewportScope | undefined {
    const scope = layer.props?.[MDV_DECK_LAYER_VIEWPORT_SCOPE];
    if (scope === "chart-shared" || scope === "per-viewport") {
        return scope;
    }
    return undefined;
}

/** Fallback for layers not yet tagged (selection id prefix, explicit viewId, viv detail id). */
export function inferDeckLayerViewportScope(layer: DeckLayerScopeInput): DeckLayerViewportScope {
    const explicit = getDeckLayerViewportScope(layer);
    if (explicit) return explicit;
    if (layer.id.startsWith("selection_")) {
        return "chart-shared";
    }
    if (typeof layer.props?.viewId === "string") {
        return "per-viewport";
    }
    return "chart-shared";
}

type CloneableDeckLayer = {
    clone: (props: Record<string, unknown>) => unknown;
};

export function tagDeckLayerViewportScope<L extends CloneableDeckLayer>(
    layer: L,
    scope: DeckLayerViewportScope,
    extraProps?: Record<string, unknown>,
): L {
    return layer.clone({
        ...extraProps,
        [MDV_DECK_LAYER_VIEWPORT_SCOPE]: scope,
    }) as L;
}

export type DeckLayerViewportFilterOptions = {
    /** Viv id suffix for the overlay detail viewport, e.g. `-#my-chartdetail-react#`. */
    overlayDetailVivId?: string;
};

/**
 * Whether `layer` should draw in `viewportId` (MDVivViewer detail view or density-grid cell).
 */
export function shouldDrawDeckLayerInViewport(
    layer: DeckLayerScopeInput,
    viewportId: string,
    options?: DeckLayerViewportFilterOptions,
): boolean {
    const scope = inferDeckLayerViewportScope(layer);
    const overlayVivId = options?.overlayDetailVivId;

    if (scope === "chart-shared") {
        if (overlayVivId && layer.id.includes(overlayVivId)) {
            return true;
        }
        return isDensityGridViewport(viewportId);
    }

    if (!isDensityGridViewport(viewportId)) {
        return overlayVivId ? layer.id.includes(overlayVivId) : false;
    }

    const layerViewId = layer.props?.viewId;
    if (typeof layerViewId === "string") {
        return layerViewId === viewportId || matchesDensityGridView(layer.id, viewportId);
    }
    return matchesDensityGridView(layer.id, viewportId);
}
