import { type RenderBackend, renderBackendSetting } from "@/lib/browserSettings";
import type { CreateDeviceProps } from "@luma.gl/core";
import { useCallback, useEffect, useState } from "react";
// Side effect: when "WebGPU call tracing" is on, wrap the WebGPU API before any device exists.
import "./webgpu_trace";

/**
 * Which deck.gl device backend the SpatialData chart renders with. WebGPU is the
 * default: MDV's own scatter layers support it (see `MdvScatterplotLayer`), and
 * contours are drawn by @intraspatial/deck's `DensityLayer`. SpatialData.js drops Viv
 * image extensions on WebGPU, so the image channel contrast/brightness (tone)
 * controls from `VivContrastExtension` have no effect there.
 *
 * A chart falls back to WebGL, for the rest of its session, when WebGPU is unavailable
 * (the probe below), when deck cannot create the device or loses it, and, if the chart
 * follows the default, when a layer fails or the canvas throws on WebGPU (see
 * `SpatialDataMDVReactComponent`). A chart set to WebGPU explicitly keeps it for those.
 *
 * The default is the "Default rendering backend" browser setting (Debug dialog; see
 * `renderBackendSetting` in `src/lib/browserSettings.ts`); a chart's own `render_backend`
 * setting overrides it, and `?webgpu` / `?webgpu=0` in the URL overrides the default for
 * one page load. The setting is observable, so charts following the default switch as
 * soon as it changes.
 *
 * deck chooses its device once, when it is created, so a change here remounts the canvas.
 */
export type DeckBackend = RenderBackend;

function urlOverride(): DeckBackend | undefined {
    if (typeof window === "undefined") return undefined;
    try {
        const param = new URLSearchParams(window.location.search).get("webgpu");
        if (param === null) return undefined;
        return param === "0" || param === "false" ? "webgl" : "webgpu";
    } catch {
        return undefined;
    }
}

const URL_OVERRIDE = urlOverride();

/** The backend for charts that leave it at the default. Observable (reads the browser setting). */
export function defaultDeckBackend(): DeckBackend {
    return URL_OVERRIDE ?? renderBackendSetting.get();
}

export function resolveDeckBackend(configured: DeckBackend | null | undefined): DeckBackend {
    return configured ?? defaultDeckBackend();
}

/** Cheap synchronous check, for deciding whether to offer WebGPU at all. */
export function browserHasWebGPU(): boolean {
    return typeof navigator !== "undefined" && "gpu" in navigator;
}

export type DeckDevice =
    | { status: "pending" }
    | {
          status: "ready";
          backend: DeckBackend;
          /** Pass as `deckProps.deviceProps`; undefined leaves deck on its WebGL default. */
          deviceProps?: CreateDeviceProps;
          /** Set when WebGPU was requested but could not be used. */
          fallbackReason?: string;
      };

const WEBGL = { status: "ready", backend: "webgl" } as const satisfies DeckDevice;

// Module-level, so every chart shares one probe and one stable `deviceProps` object.
let webgpuDevice: Promise<DeckDevice> | null = null;

/** Every limit the adapter reports, at its own value: what luma's "max" level requests. */
function adapterLimits(limits: GPUSupportedLimits): Record<string, number> {
    const required: Record<string, number> = {};
    // The limits are getters on the prototype, so `for...in` rather than `Object.entries`.
    for (const name in limits) {
        const value: unknown = Reflect.get(limits, name);
        if (typeof value === "number") required[name] = value;
    }
    return required;
}

async function probeWebGPU(): Promise<DeckDevice> {
    const fallback = (fallbackReason: string): DeckDevice => ({ ...WEBGL, fallbackReason });
    if (!browserHasWebGPU()) return fallback("this browser does not support WebGPU");
    try {
        const adapter = await navigator.gpu.requestAdapter();
        if (!adapter) return fallback("no WebGPU adapter is available");
        // luma's default "core" level gets WebGPU's spec-minimum limits, so a 10M-point
        // scatter (480MB of instance attributes) overflows the 256MB `maxBufferSize`.
        // "max" asks for every limit and feature at the adapter's own value. That also
        // covers DensityLayer's faster path, which blends into r32float
        // ("float32-blendable"); without it (Firefox, for one) the layer falls back to an
        // exact compute splat.
        // An adapter can still refuse a device. Ask for the limits deck will ask for (the
        // features are the adapter's own, so cannot be refused), so that failure falls
        // back here instead of leaving deck with no device.
        const probe = await adapter.requestDevice({ requiredLimits: adapterLimits(adapter.limits) });
        probe.destroy();
        const { webgpuAdapter } = await import("@luma.gl/webgpu");
        return {
            status: "ready",
            backend: "webgpu",
            deviceProps: {
                type: "webgpu",
                adapters: [webgpuAdapter],
                featureLevel: "max",
            },
        };
    } catch (error) {
        return fallback(error instanceof Error ? error.message : String(error));
    }
}

/**
 * Resolve the requested backend to deck `deviceProps`, falling back to WebGL.
 *
 * `fallBack` is for a WebGPU device deck could not create although the probe passed
 * (luma's device setup, or the canvas's WebGPU context): it puts this chart on WebGL
 * until its backend changes. It is per chart because a canvas context can fail alone.
 */
export function useDeckDevice(backend: DeckBackend): { device: DeckDevice; fallBack: (reason: string) => void } {
    const [webgpu, setWebgpu] = useState<DeckDevice>({ status: "pending" });
    const [failure, setFailure] = useState<string | null>(null);
    // A backend change gives WebGPU another try (React's reset-on-prop-change pattern).
    const [lastBackend, setLastBackend] = useState(backend);
    if (backend !== lastBackend) {
        setLastBackend(backend);
        setFailure(null);
    }
    useEffect(() => {
        if (backend !== "webgpu") return;
        let cancelled = false;
        webgpuDevice ??= probeWebGPU();
        webgpuDevice.then((device) => {
            if (!cancelled) setWebgpu(device);
        });
        return () => {
            cancelled = true;
        };
    }, [backend]);
    const fallBack = useCallback((reason: string) => setFailure(reason), []);
    if (backend !== "webgpu") return { device: WEBGL, fallBack };
    if (failure) return { device: { ...WEBGL, fallbackReason: failure }, fallBack };
    return { device: webgpu, fallBack };
}
