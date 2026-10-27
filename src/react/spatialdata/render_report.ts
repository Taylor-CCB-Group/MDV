import type { DeviceInfo } from "@luma.gl/core";
import type { CanvasError, LayerFailure } from "./layer_failures";

/**
 * A SpatialData chart's render problems as plain text, for the status indicator's "Copy
 * details". Plain text rather than JSON (as `DebugErrorComponent` copies), so stacks keep
 * their line breaks when pasted into an issue.
 */

/** One line describing the GPU, from luma's device info. */
export function describeGpu(info: DeviceInfo): string {
    const details = [info.gpu, info.gpuType, info.gpuArchitecture, info.gpuBackend, info.featureLevel]
        .filter((part) => part && part !== "unknown")
        .join(", ");
    const name = [info.vendor, info.renderer].filter(Boolean).join(" ");
    return `${info.type}: ${name || "unknown GPU"}${details ? ` (${details})` : ""}`;
}

export type RenderReportInput = {
    renderer: string;
    gpu: string | null;
    userAgent: string;
    url: string;
    failures: Array<{ name: string; failure: LayerFailure }>;
    canvasErrors: readonly CanvasError[];
    /** `webgpuTraceReport()`: null while tracing is off. */
    trace?: string | null;
};

function indent(text: string) {
    return text
        .trimEnd()
        .split("\n")
        .map((line) => `    ${line}`)
        .join("\n");
}

export function formatRenderReport(input: RenderReportInput): string {
    const lines = [
        "MDV SpatialData chart: rendering problems",
        `Renderer: ${input.renderer}`,
        `GPU: ${input.gpu ?? "not known yet"}`,
        `Browser: ${input.userAgent}`,
        `Page: ${input.url}`,
    ];
    for (const { name, failure } of input.failures) {
        lines.push("", `Layer turned off: ${name} (${failure.source}, ×${failure.count})`, indent(failure.reason));
        if (failure.detail) lines.push("  Reported as:", indent(failure.detail));
        if (failure.stack) lines.push("  Stack:", indent(failure.stack));
    }
    for (const error of input.canvasErrors) {
        lines.push("", `Rendering error (×${error.count}):`, indent(error.message));
        if (error.stack) lines.push("  Stack:", indent(error.stack));
    }
    if (input.failures.length === 0 && input.canvasErrors.length === 0) lines.push("", "No rendering problems recorded.");
    lines.push("", input.trace ? `WebGPU call trace:\n${indent(input.trace)}` : "WebGPU call tracing: off (Debug dialog > Browser settings)");
    return `${lines.join("\n")}\n`;
}
