import type { CreateDeviceProps } from "@luma.gl/core";
import { useEffect, useState } from "react";

/**
 * Which deck.gl device backend the SpatialData chart renders with. WebGPU is
 * experimental: MDV's own scatter layers support it (see `MdvScatterplotLayer`),
 * but on WebGPU contours lose their contour lines, and SpatialData.js drops Viv
 * image extensions, so the image channel contrast/brightness (tone) controls from
 * `VivContrastExtension` have no effect.
 *
 * The default comes from a dev flag; a chart's `render_backend` setting overrides it.
 *
 *   localStorage.MDV_WEBGPU = "1"   // then reload; sticky, like MDV_SPATIAL_PERF
 *   delete localStorage.MDV_WEBGPU  // then reload
 *   ?webgpu  /  ?webgpu=0           // per page load, overrides localStorage
 *
 * deck chooses its device once, when it is created, so a change here remounts the canvas.
 */
export type DeckBackend = "webgl" | "webgpu";

const FLAG = "MDV_WEBGPU";

function readFlag(): boolean {
    if (typeof window === "undefined") return false;
    try {
        const param = new URLSearchParams(window.location.search).get("webgpu");
        if (param !== null) return param !== "0" && param !== "false";
        return window.localStorage?.getItem(FLAG) === "1";
    } catch {
        return false;
    }
}

const FLAG_SET = readFlag();
if (FLAG_SET) {
    console.warn(
        `[MDV] WebGPU is the default backend for SpatialData charts (${FLAG} flag or ?webgpu). ` +
            `This is experimental; \`delete localStorage.${FLAG}\` and reload to turn it off.`,
    );
}

export function resolveDeckBackend(configured: DeckBackend | null | undefined): DeckBackend {
    return configured ?? (FLAG_SET ? "webgpu" : "webgl");
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
        return { status: "ready", backend: "webgpu", deviceProps: { type: "webgpu", adapters: [webgpuAdapter] } };
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
