import type { Layer } from "@deck.gl/core";
import { useMemo } from "react";
import {
    cloneDeckLayerForRender,
    cloneLayerForChartArrayGrid,
} from "./densityGridUtils";
import {
    buildChartArrayDeckLayers,
    type BuildChartArrayDeckLayersOptions,
} from "./chartArrayDeckLayers";

type CloneableDeckLayer = Layer & {
    clone: (props: Record<string, unknown>) => Layer;
};

const templateInstanceIds = new WeakMap<Layer, number>();
let nextTemplateInstanceId = 0;

function getTemplateInstanceId(layer: Layer): number {
    const existing = templateInstanceIds.get(layer);
    if (existing !== undefined) return existing;
    const id = nextTemplateInstanceId++;
    templateInstanceIds.set(layer, id);
    return id;
}

/** Fingerprint template layer props that affect deck initialization. */
export function getDeckLayerTemplateSignature(layer: Layer): string {
    const data = layer.props.data;
    let dataKey = "none";
    if (typeof data === "object" && data !== null && "length" in data) {
        dataKey = `len:${String((data as { length: number }).length)}`;
    } else if (Array.isArray(data)) {
        dataKey = `arr:${data.length}`;
    } else if (data instanceof Uint32Array || data instanceof Float32Array) {
        dataKey = `typed:${data.length}`;
    }
    return `${layer.id}|vis=${String(layer.props.visible)}|${dataKey}|inst=${getTemplateInstanceId(layer)}`;
}

export function cloneDeckTemplateLayer(
    template: Layer,
    mode: "overlay" | "chart-array-grid",
    extraProps: Record<string, unknown> = {},
): Layer {
    const cloneable = template as CloneableDeckLayer;
    if (mode === "chart-array-grid") {
        return cloneLayerForChartArrayGrid(cloneable, extraProps);
    }
    return cloneDeckLayerForRender(cloneable, extraProps);
}

type ClonedDeckLayersOptions = {
    mode: "overlay" | "chart-array-grid";
    extraPropsForTemplate?: (template: Layer) => Record<string, unknown>;
};

/**
 * Clone spatial templates for deck render. Templates from useSpatialLayers must
 * never be passed to deck — each mount gets fresh clone instances (new JS objects,
 * stable ids) so deck.gl does not assert on finalized layer reuse.
 */
export function useClonedDeckLayersForDeck(
    templates: readonly (Layer | null | undefined)[],
    deckMountEpoch: number,
    options: ClonedDeckLayersOptions,
    active = true,
): Layer[] {
    const signatures = templates.map((layer) =>
        layer ? getDeckLayerTemplateSignature(layer) : "",
    );

    // biome-ignore lint/correctness/useExhaustiveDependencies: signatures + deckMountEpoch capture template churn; omit templates/extraProps (fresh clones each rebuild)
    return useMemo(() => {
        if (!active) {
            return [];
        }
        return templates
            .filter((layer): layer is Layer => layer != null)
            .map((template) => {
                const extraProps = options.extraPropsForTemplate?.(template) ?? {};
                return cloneDeckTemplateLayer(template, options.mode, extraProps);
            });
    }, [active, deckMountEpoch, options.mode, ...signatures]);
}

export function buildDeckLayerCacheKey(parts: readonly string[]): string {
    return parts.filter(Boolean).join("\u0001");
}

/** Chart-array grid deck layers — rebuilt when cacheKey or deckMountEpoch changes. */
export function useBuiltChartArrayDeckLayers(
    options: BuildChartArrayDeckLayersOptions,
    cacheKey: string,
    deckMountEpoch: number,
): Layer[] {
    // biome-ignore lint/correctness/useExhaustiveDependencies: cacheKey + deckMountEpoch intentionally subsume options churn
    return useMemo(() => buildChartArrayDeckLayers(options), [cacheKey, deckMountEpoch]);
}

/** @deprecated Use {@link useClonedDeckLayersForDeck} */
export const useStableClonedDeckLayers = useClonedDeckLayersForDeck;

/** @deprecated Use {@link useBuiltChartArrayDeckLayers} */
export const useStableBuiltChartArrayDeckLayers = useBuiltChartArrayDeckLayers;
