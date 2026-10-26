import { type RenderBackend, renderBackendSetting } from "@/lib/browserSettings";
import type { CreateDeviceProps } from "@luma.gl/core";
import { useEffect, useState } from "react";

/**
 * Which deck.gl device backend the SpatialData chart renders with. WebGPU is
 * experimental: MDV's own scatter layers support it (see `MdvScatterplotLayer`), and
 * contours are drawn by @intraspatial/deck's `DensityLayer`. SpatialData.js drops Viv
 * image extensions on WebGPU, so the image channel contrast/brightness (tone)
 * controls from `VivContrastExtension` have no effect there.
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
if ((URL_OVERRIDE ?? renderBackendSetting.get()) === "webgpu") {
    console.warn(
        "[MDV] WebGPU is the default backend for SpatialData charts in this browser. This is " +
            "experimental; change 'Default rendering backend' in the Debug dialog to turn it off.",
    );
}

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

async function probeWebGPU(): Promise<DeckDevice> {
    const fallback = (fallbackReason: string): DeckDevice => ({ ...WEBGL, fallbackReason });
    if (!browserHasWebGPU()) return fallback("this browser does not support WebGPU");
    try {
        const adapter = await navigator.gpu.requestAdapter();
        if (!adapter) return fallback("no WebGPU adapter is available");
        const { webgpuAdapter } = await import("@luma.gl/webgpu");
        return {
            status: "ready",
            backend: "webgpu",
            deviceProps: {
                type: "webgpu",
                adapters: [webgpuAdapter],
                // DensityLayer's faster path blends into r32float. Optional: without it
                // (Firefox, for one) the layer falls back to an exact compute splat.
                optionalFeatures: ["float32-blendable"],
            },
        };
    } catch (error) {
        return fallback(error instanceof Error ? error.message : String(error));
    }
}

/** Resolve the requested backend to deck `deviceProps`, falling back to WebGL. */
export function useDeckDevice(backend: DeckBackend): DeckDevice {
    const [webgpu, setWebgpu] = useState<DeckDevice>({ status: "pending" });
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
    return backend === "webgpu" ? webgpu : WEBGL;
}
