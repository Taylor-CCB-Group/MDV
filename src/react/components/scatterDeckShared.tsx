import type { Layer } from "@deck.gl/core";
import { action } from "mobx";
import { type ReactNode, useCallback, useEffect, useMemo, useRef } from "react";
import { ErrorBoundary } from "react-error-boundary";
import type { DataColumn, NumberDataType } from "@/charts/charts";
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
import type { DeckScatterConfig, DeckScatterReact } from "./DeckScatterReactWrapper";
import ErrorComponentReactWrapper from "./ErrorComponentReactWrapper";

/**
 * What the deck scatter plot and its density grid share: fitting the view to the filtered
 * rows, and the device deck draws with.
 */

export const errorMessage = (error: unknown) => (error instanceof Error ? error.message : String(error));

//todo this should be in a common place etc.
const colMid = ({ minMax }: DataColumn<NumberDataType>) => minMax[0] + (minMax[1] - minMax[0]) / 2;

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
        //there is also cx.minMax, cy.minMax - but not for filtered indices
        let minX = Number.POSITIVE_INFINITY;
        let maxX = Number.NEGATIVE_INFINITY;
        let minY = Number.POSITIVE_INFINITY;
        let maxY = Number.NEGATIVE_INFINITY;
        // todo z... probably refactor for less repetition
        // ^ maybe we could have properties on the columns that are the min/max of the filtered data
        for (let i = 0; i < data.length; i++) {
            const d = data[i];
            const x = cx.data[d];
            const y = cy.data[d];
            if (!Number.isFinite(x) || !Number.isFinite(y)) {
                console.warn("undefined data in scatterplot");
                continue;
            }
            if (x < minX) minX = x;
            if (x > maxX) maxX = x;
            if (y < minY) minY = y;
            if (y > maxY) maxY = y;
        }
        action(() => {
            chart.ignoreStateUpdate = true;
            // the 3D branch can return or throw early: never leave updates suppressed
            try {
                chart.pendingRecenter = false;

                const scaling = config.dimension === "2d" ? (config.axis_scaling ?? "equal") : "equal";
                config.viewState = fitOrthographicView({ minX, maxX, minY, maxY }, width, height, scaling);
                // the stored view (mobx's copy), to tell later whether anyone has moved it
                lastFit.current = { width, height, viewState: config.viewState };
                // any kind of interpolator would probably mean changing how we handle viewState
                // transitionInterpolator: new FlyToInterpolator({speed: 1.2}),
                if (config.dimension === "3d") {
                    const vs = config.viewState;
                    vs.rotationOrbit = vs.rotationOrbit || 0;
                    vs.rotationX = vs.rotationX || 0;
                    if (!config.zoom_on_filter) {
                        if (!cz) {
                            console.warn("3D viewState requires z coordinate");
                            return;
                        }
                        const params = [cx, cy, cz];
                        if (!allNumeric(params)) {
                            throw new Error("3D viewState requires numeric coordinates");
                        }
                        vs.target = params.map(colMid) as [number, number, number];
                    }
                }
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
