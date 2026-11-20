import type { Layer } from "@deck.gl/core";
import { action } from "mobx";
import { type ReactNode, useCallback, useEffect, useMemo, useRef } from "react";
import { ErrorBoundary } from "react-error-boundary";
import { allNumeric } from "@/lib/columnTypeHelpers";
import {
    type DeckDevice,
    type DeckErrorHandlers,
    resolveDeckBackend,
    useDeckDevice,
    useDeckDeviceProps,
} from "@/react/spatialdata/deck_backend";
import { LayerFailureTracker } from "@/react/spatialdata/layer_failures";
import { fitOrthographicView } from "../axisScaling";
import { useChart } from "../context";
import { useConfig, useParamColumns } from "../hooks";
import { fitOrbitView, validAxisScales } from "../orbitView";
import type { DeckScatterConfig, DeckScatterReact } from "./DeckScatterReactWrapper";
import ErrorComponentReactWrapper from "./ErrorComponentReactWrapper";

/**
 * What the deck scatter plot and its density grid share: fitting the view to the filtered
 * rows, and the device deck draws with.
 */

export const errorMessage = (error: unknown) => (error instanceof Error ? error.message : String(error));

/** The smallest and largest finite value of each column among `rows`; undefined if there are none. */
function filteredBounds(rows: Uint32Array, columns: readonly ArrayLike<number>[]) {
    const n = columns.length;
    const min = new Array<number>(n).fill(Number.POSITIVE_INFINITY);
    const max = new Array<number>(n).fill(Number.NEGATIVE_INFINITY);
    let skipped = 0;
    nextRow: for (let r = 0; r < rows.length; r++) {
        const row = rows[r];
        for (let c = 0; c < n; c++) {
            if (!Number.isFinite(columns[c][row])) {
                skipped++;
                continue nextRow;
            }
        }
        for (let c = 0; c < n; c++) {
            const v = columns[c][row];
            if (v < min[c]) min[c] = v;
            if (v > max[c]) max[c] = v;
        }
    }
    if (skipped) console.warn(`${skipped} rows with missing coordinates in scatterplot`);
    if (!Number.isFinite(min[0])) return undefined;
    return { min, max };
}

// const margin = { top: 10, right: 10, bottom: 40, left: 60 };
/** todo this should be common for viv / scatter_state, pending refactor
 * there should be hooks getting the range of a filtered column
 * so that multiple charts can re-use the computation (but if not used, it's not computed)
 */
export function useZoomOnFilter(
    data: Uint32Array,
    width: number,
    height: number,
    {
        refitOnResize = false,
    }: {
        /**
         * Fit again when the size changes, while the view is still the last fit (nobody has
         * panned or zoomed since). For a density grid, whose cells are first measured before
         * the grid's layout settles.
         */
        refitOnResize?: boolean;
    } = {},
) {
    const [cx, cy, cz] = useParamColumns();
    // const data = useFilteredIndices(); //<< does calling this multiple times mean wasting space?
    const config = useConfig<DeckScatterConfig>();
    const chart = useChart() as any;
    const { pendingRecenter } = chart;
    const lastFit = useRef<{ width: number; height: number; viewState: unknown } | null>(null);
    useEffect(() => {
        if (chart.ignoreStateUpdate) return;
        if (data.length === 0) return; // [0, 0, 1, 1];
        if (!(width > 0 && height > 0)) return;
        const fit = lastFit.current;
        const resized =
            refitOnResize &&
            fit !== null &&
            fit.viewState === config.viewState &&
            (fit.width !== width || fit.height !== height);
        if (!pendingRecenter && !config.zoom_on_filter && !resized) return;
        if (config.dimension === "3d") {
            if (!cz || !allNumeric([cx, cy, cz])) {
                console.warn("3D viewState requires numeric x, y and z columns");
                return;
            }
            const scales = validAxisScales(config.axis_scales);
            const bounds = filteredBounds(data, [cx.data, cy.data, cz.data]);
            if (!bounds) return;
            // zoom on filter frames the rows where they are, from where the camera is
            const previous = config.viewState;
            const rotation = pendingRecenter
                ? undefined
                : { rotationOrbit: previous.rotationOrbit ?? 0, rotationX: previous.rotationX ?? 0 };
            action(() => {
                chart.ignoreStateUpdate = true;
                try {
                    chart.pendingRecenter = false;
                    config.viewState = fitOrbitView(
                        {
                            min: bounds.min.map((v, i) => v * scales[i]),
                            max: bounds.max.map((v, i) => v * scales[i]),
                        },
                        width,
                        height,
                        rotation,
                    );
                    lastFit.current = { width, height, viewState: config.viewState };
                } finally {
                    chart.ignoreStateUpdate = false;
                }
            })();
            return;
        }
        const bounds = filteredBounds(data, [cx.data, cy.data]);
        if (!bounds) return;
        const [minX, minY] = bounds.min;
        const [maxX, maxY] = bounds.max;
        action(() => {
            chart.ignoreStateUpdate = true;
            try {
                chart.pendingRecenter = false;
                const scaling = config.dimension === "2d" ? (config.axis_scaling ?? "equal") : "equal";
                config.viewState = fitOrthographicView({ minX, maxX, minY, maxY }, width, height, scaling);
                // the stored view (mobx's copy), to tell later whether anyone has moved it
                lastFit.current = { width, height, viewState: config.viewState };
                // any kind of interpolator would probably mean changing how we handle viewState
                // transitionInterpolator: new FlyToInterpolator({speed: 1.2}),
            } finally {
                chart.ignoreStateUpdate = false;
            }
        })();
    }, [data, cx, cy, cz, width, height, config, pendingRecenter, chart, refitOnResize]);
}

/**
 * The device deck draws with (see `deck_backend.ts`), and what to do when it fails. On
 * WebGPU, a chart following the default backend moves to WebGL when a layer, the device or
 * the canvas fails; one set to WebGPU explicitly stays on it, with the errors in the console
 * (each kind once: a WebGPU validation error otherwise repeats every frame).
 */
export function useScatterDeckDevice(layerIds: readonly string[]) {
    const chart = useChart<DeckScatterConfig, DeckScatterReact>();
    const config = useConfig<DeckScatterConfig>();
    const { device, fallBack } = useDeckDevice(resolveDeckBackend(config.render_backend));
    const activeBackend = device.status === "ready" ? device.backend : null;
    const fallbackReason = device.status === "ready" ? device.fallbackReason : undefined;
    useEffect(() => {
        chart.setActiveBackend(activeBackend, fallbackReason);
        return () => chart.setActiveBackend(null);
    }, [chart, activeBackend, fallbackReason]);
    const layerIdsRef = useRef(layerIds);
    layerIdsRef.current = layerIds;
    /** Fall back to WebGL, if this chart is on WebGPU because it follows the default. */
    const failedOnWebGPU = useCallback(
        (reason: string) => {
            // after the switch, the old WebGPU device's errors are stale
            if (chart.activeBackend !== "webgpu" || chart.config.render_backend != null) return;
            console.error(`[MDV] ${reason}; using WebGL`);
            fallBack(reason);
        },
        [chart, fallBack],
    );
    const failureTracker = useMemo(
        () =>
            new LayerFailureTracker({
                entryIds: () => layerIdsRef.current,
                fail: (layerId, failure) => failedOnWebGPU(`${layerId} failed on WebGPU: ${failure.reason}`),
                canvas: (message) => failedOnWebGPU(`the canvas failed on WebGPU: ${message}`),
            }),
        [failedOnWebGPU],
    );
    useEffect(() => () => failureTracker.dispose(), [failureTracker]);
    const handlers = useMemo<DeckErrorHandlers>(
        () => ({
            device: (error) => failureTracker.device(error),
            layer: (error, layerId) => failureTracker.layer(error, layerId),
        }),
        [failureTracker],
    );
    const deckProps = useDeckDeviceProps(device, fallBack, handlers);
    return { device, failedOnWebGPU, ...deckProps };
}

/**
 * A deck layer instance belongs to the first Deck that draws it: its models and buffers are
 * on that Deck's device. These layers are memoised outside the canvas, so a Deck remounted
 * for another backend would inherit the old device's state and fail to draw. Give each
 * device its own copies, one per layer instance.
 */
export function useLayersForDevice(layers: Layer[], device: DeckDevice): Layer[] {
    // biome-ignore lint/correctness/useExhaustiveDependencies: a new device is the trigger
    const copies = useMemo(() => new WeakMap<Layer, Layer>(), [device]);
    return layers.map((layer) => {
        const cached = copies.get(layer);
        if (cached) return cached;
        const copy: Layer = layer.clone({});
        copies.set(layer, copy);
        return copy;
    });
}

/**
 * Shows the deck canvas once its device is ready. Deck picks its device once, when created,
 * so the canvas remounts on a backend change; a canvas that fails on WebGPU falls back.
 */
export function ScatterDeckBoundary({
    device,
    failedOnWebGPU,
    children,
}: {
    device: DeckDevice;
    failedOnWebGPU: (reason: string) => void;
    children: ReactNode;
}) {
    if (device.status === "pending") return <div className="h-full w-full p-2">Starting WebGPU…</div>;
    return (
        <ErrorBoundary
            key={device.backend}
            onError={(error) => {
                if (device.backend === "webgpu") failedOnWebGPU(`the canvas failed on WebGPU: ${errorMessage(error)}`);
            }}
            fallbackRender={({ error }) => (
                <ErrorComponentReactWrapper
                    error={{
                        message: errorMessage(error),
                        ...(error instanceof Error && error.stack ? { stack: error.stack } : {}),
                    }}
                    title="Error showing scatter plot canvas"
                />
            )}
        >
            {children}
        </ErrorBoundary>
    );
}
