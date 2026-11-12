import type { DeckGLRef, Layer, LayersList } from "deck.gl";
import { type RefObject, useCallback, useRef, useState } from "react";
import { recordSpatialLoaded } from "./perf";

/**
 * Finds when the SpatialData.js canvas has drawn every layer fully loaded, and reports
 * how long that took through `recordSpatialLoaded` in `perf.ts`.
 *
 * deck.gl's `layer.isLoaded` is only true once a layer's data has arrived, and for the
 * tile layers that draw images, once every tile in view has loaded. So the first frame
 * drawn while every layer reports loaded is when the user sees the finished picture.
 *
 * - `first`: from when the chart opened (`useMarkSpatialChartOpened`) to that frame.
 * - `settle`: from the last frame whose viewport moved to that frame. Only recorded when
 *   the move left some layer loading, so moving over tiles already in memory adds nothing.
 *
 * The viewport is read from deck itself after each frame, because the chart's view
 * state prop does not change on every pan or zoom.
 */

const openedAt = new WeakMap<object, number>();

/** Call in the chart's root component, so `first` includes opening the store. */
export function useMarkSpatialChartOpened(chart: object) {
    useState(() => {
        if (!openedAt.has(chart)) openedAt.set(chart, performance.now());
        return null;
    });
}

function flattenLayers(layers: LayersList, out: Layer[] = []): Layer[] {
    for (const layer of layers) {
        if (Array.isArray(layer)) flattenLayers(layer, out);
        else if (layer) out.push(layer as Layer);
    }
    return out;
}

type Pending = { kind: "first" | "settle"; since: number; waited: boolean };

export function useSpatialLoadTiming(
    deckRef: RefObject<DeckGLRef | null>,
    chart: object,
    describe: () => Record<string, string>,
) {
    const pending = useRef<Pending | null>({
        kind: "first",
        since: openedAt.get(chart) ?? performance.now(),
        waited: true,
    });
    const lastViewport = useRef<string | null>(null);

    // Runs after every frame, and checks even while nothing records, so a timing that
    // starts recording later does not report a load that finished long before. The
    // work is a viewport comparison, plus a few layers while something is loading.
    return useCallback(() => {
        const deck = deckRef.current?.deck;
        if (!deck) return;
        const now = performance.now();

        const viewport = deck.getViewports()[0];
        if (viewport) {
            const key = `${viewport.width}x${viewport.height} ${viewport.zoom} ${viewport.center}`;
            const moved = lastViewport.current !== null && key !== lastViewport.current;
            lastViewport.current = key;
            // Opening the chart frames the view; that wait belongs to `first`.
            if (moved && pending.current?.kind !== "first") {
                pending.current = { kind: "settle", since: now, waited: false };
            }
        }

        const current = pending.current;
        if (!current) return;
        const drawn = flattenLayers(deck.props.layers);
        if (drawn.length === 0) return;
        if (!drawn.every((layer) => layer.isLoaded)) {
            current.waited = true;
            return;
        }
        pending.current = null;
        if (current.waited) recordSpatialLoaded(current.kind, now - current.since, describe());
    }, [deckRef, describe]);
}
