import type { DeckDevice } from "@/react/spatialdata/deck_backend";
import { renderStackEntryDisplayName } from "@/react/spatialdata/render_stack_display";
import { Popover } from "@mui/material";
import type { SpatialData } from "@spatialdata/core";
import { observer } from "mobx-react-lite";
import { type ReactNode, useState } from "react";
import { createPortal } from "react-dom";
import { useChart } from "../context";
import CustomTooltip from "./CustomTooltip";
import { RasterIssueList, rasterEntriesKey, useRasterIssues } from "./RasterIssues";
import type { SpatialDataMdvReact, SpatialDataMdvReactConfig } from "./SpatialDataMDVReact";

type ReadyDevice = Extract<DeckDevice, { status: "ready" }>;

function backendLabel(device: ReadyDevice) {
    if (device.backend === "webgpu") return "WebGPU";
    if (device.fallbackReason) return `WebGL (WebGPU unavailable: ${device.fallbackReason})`;
    return "WebGL";
}

function Section({ title, children }: { title: string; children: ReactNode }) {
    return (
        <section className="flex flex-col gap-1">
            <h3 className="text-xs font-semibold">{title}</h3>
            {children}
        </section>
    );
}

/**
 * The chart's status, in its title bar rather than over the canvas: a red warning with a count
 * while there are render errors nobody has dismissed, a plain notice while layers stay off or
 * image data is slow, and a popover with the details — failed layers, unattributed render
 * errors, slow image data, and which renderer is in use.
 */
const SpatialStatusIndicator = observer(function SpatialStatusIndicator({
    spatialData,
    device,
}: {
    spatialData: SpatialData;
    device: ReadyDevice;
}) {
    const chart = useChart<SpatialDataMdvReactConfig, SpatialDataMdvReact>();
    const [anchor, setAnchor] = useState<HTMLElement | null>(null);
    const rasterIssues = useRasterIssues(spatialData, rasterEntriesKey(chart.config.renderStack));
    const entries = chart.config.renderStack?.entries ?? [];
    const failures = [...chart.layerFailures].map(([entryId, failure]) => {
        const entry = entries.find((item) => item.id === entryId);
        return { entryId, failure, name: entry ? renderStackEntryDisplayName(entry) : entryId };
    });
    const errorCount = failures.filter((item) => !item.failure.dismissed).length + chart.canvasErrors.length;
    const hasNotice = failures.length > 0 || rasterIssues.length > 0;
    const spatialDataUrl = spatialData.url;

    const summary = [
        errorCount > 0 && `${errorCount} rendering ${errorCount === 1 ? "problem" : "problems"}`,
        errorCount === 0 &&
            failures.length > 0 &&
            `${failures.length} ${failures.length === 1 ? "layer" : "layers"} turned off`,
        rasterIssues.length > 0 && "slow image data",
        `rendering with ${backendLabel(device)}`,
    ]
        .filter(Boolean)
        .join(" · ");

    return createPortal(
        <>
            <CustomTooltip tooltipText={summary}>
                <button
                    type="button"
                    aria-label={`Chart status: ${summary}`}
                    className="flex cursor-pointer items-center gap-0.5 border-0 bg-transparent p-0 text-inherit"
                    style={errorCount > 0 ? { color: "var(--icon_color_error)" } : undefined}
                    onClick={(event) => setAnchor(event.currentTarget)}
                >
                    <i
                        className={`ciview-chart-icon fas ${
                            errorCount > 0
                                ? "fa-exclamation-triangle"
                                : hasNotice
                                  ? "fa-exclamation-circle"
                                  : "fa-check-circle"
                        }`}
                        style={errorCount === 0 && !hasNotice ? { opacity: 0.6 } : undefined}
                    />
                    {errorCount > 0 && <span className="text-xs">{errorCount}</span>}
                </button>
            </CustomTooltip>
            <Popover
                open={anchor !== null}
                anchorEl={anchor}
                onClose={() => setAnchor(null)}
                anchorOrigin={{ vertical: "bottom", horizontal: "right" }}
                transformOrigin={{ vertical: "top", horizontal: "right" }}
            >
                <div className="flex max-h-[60vh] max-w-[26rem] flex-col gap-3 overflow-y-auto p-3 text-xs">
                    {failures.length > 0 && (
                        <Section title="Layers turned off after failing">
                            {failures.map(({ entryId, failure, name }) => (
                                <p
                                    key={entryId}
                                    className="line-clamp-3 break-words"
                                    title={[failure.reason, failure.detail].filter(Boolean).join("\n\n")}
                                >
                                    <span className="font-medium">{name}</span>: {failure.reason}
                                </p>
                            ))}
                            <button
                                type="button"
                                className="cursor-pointer self-start border-0 bg-transparent p-0 text-inherit underline"
                                onClick={() => {
                                    setAnchor(null);
                                    chart.openLayerDialog();
                                }}
                            >
                                Open the layer dialog
                            </button>
                        </Section>
                    )}
                    {chart.canvasErrors.length > 0 && (
                        <Section title="Rendering errors">
                            {chart.canvasErrors.map((error) => (
                                <p key={error.message} className="line-clamp-3 break-words" title={error.message}>
                                    {error.message}
                                    {error.count > 1 && ` (×${error.count})`}
                                </p>
                            ))}
                            <p className="opacity-70">Details are in the browser console.</p>
                        </Section>
                    )}
                    {errorCount > 0 && (
                        <button
                            type="button"
                            className="cursor-pointer self-start rounded border border-current bg-transparent px-2 py-0.5 text-inherit"
                            onClick={() => chart.dismissRenderErrors()}
                        >
                            Dismiss errors
                        </button>
                    )}
                    {rasterIssues.length > 0 && (
                        <Section title="Slow image data">
                            <RasterIssueList issues={rasterIssues} />
                        </Section>
                    )}
                    <Section title="Renderer">
                        <p>{backendLabel(device)}</p>
                    </Section>
                    {spatialDataUrl && (
                        <a
                            href={`https://taylor-ccb-group.github.io/SpatialData.js/docs/demo/?url=${encodeURIComponent(spatialDataUrl)}`}
                            target="_blank"
                            rel="noopener noreferrer"
                            className="text-inherit underline"
                        >
                            Open in SpatialData.js
                        </a>
                    )}
                </div>
            </Popover>
        </>,
        chart.statusIndicatorHost,
    );
});

export default SpatialStatusIndicator;
