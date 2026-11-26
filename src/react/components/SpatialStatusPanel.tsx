import { Button } from "@mui/material";
import { observer } from "mobx-react-lite";
import { type ReactNode, useState } from "react";
import { formatRenderReport } from "@/react/spatialdata/render_report";
import { webgpuTraceReport } from "@/react/spatialdata/webgpu_trace";
import EntryStatusDetails from "./EntryStatusDetails";
import type { SpatialDataMdvReact } from "./SpatialDataMDVReact";

function Section({ title, children }: { title: string; children: ReactNode }) {
    return (
        <section className="flex flex-col gap-1">
            <h3 className="text-sm font-semibold">{title}</h3>
            {children}
        </section>
    );
}

/**
 * The chart's status, in its debug dialog (`getDebugPanel`): failed layers, render errors
 * that could not be tied to a layer, what sd.js reports about each layer (load errors, and
 * notices such as slow image data), and which renderer is in use. The chart's menu pip
 * (`statusPip`) points here.
 */
const SpatialStatusPanel = observer(function SpatialStatusPanel({ chart }: { chart: SpatialDataMdvReact }) {
    const [copyState, setCopyState] = useState<"idle" | "copied" | "failed">("idle");
    const { failures, statuses, errorCount } = chart.renderStatus;
    const { canvasErrors, spatialDataUrl, retryEntry } = chart;
    const failedToLoad = statuses.filter(({ status }) => status.errors.length > 0);
    const withNotices = statuses.filter(({ status }) => status.notices.length > 0);

    const copyDetails = async () => {
        const report = formatRenderReport({
            renderer: chart.rendererLabel,
            gpu: chart.gpuDescription,
            userAgent: navigator.userAgent,
            url: window.location.href,
            failures: failures.map(({ name, failure }) => ({ name, failure })),
            canvasErrors,
            loadErrors: failedToLoad.map(({ name, status }) => ({ name, errors: status.errors })),
            trace: webgpuTraceReport(),
        });
        try {
            await navigator.clipboard.writeText(report);
            setCopyState("copied");
        } catch (error) {
            // e.g. no clipboard permission: the console still has it
            console.info(report);
            console.warn("Could not copy the render report:", error);
            setCopyState("failed");
        }
        setTimeout(() => setCopyState("idle"), 2000);
    };

    return (
        <div className="flex flex-col gap-3 text-sm">
            {failures.length > 0 && (
                <Section title="Layers turned off after failing">
                    {failures.map(({ entryId, failure, name }) => (
                        <p
                            key={entryId}
                            className="line-clamp-3 wrap-break-word"
                            title={[failure.reason, failure.detail].filter(Boolean).join("\n\n")}
                        >
                            <span className="font-medium">{name}</span>: {failure.reason}
                        </p>
                    ))}
                    <button
                        type="button"
                        className="cursor-pointer self-start border-0 bg-transparent p-0 text-inherit underline"
                        onClick={() => chart.openLayerDialog()}
                    >
                        Open the layer dialog
                    </button>
                </Section>
            )}
            {canvasErrors.length > 0 && (
                <Section title="Rendering errors">
                    {canvasErrors.map((error) => (
                        <p key={error.message} className="line-clamp-3 wrap-break-word" title={error.message}>
                            {error.message}
                            {error.count > 1 && ` (×${error.count})`}
                        </p>
                    ))}
                </Section>
            )}
            {(errorCount > 0 || failures.length > 0) && (
                <div className="flex flex-wrap items-center gap-2">
                    <Button variant="outlined" size="small" onClick={copyDetails}>
                        {copyState === "copied" ? "Copied" : "Copy details"}
                    </Button>
                    {errorCount > 0 && (
                        <Button variant="outlined" size="small" onClick={() => chart.dismissRenderErrors()}>
                            Dismiss errors
                        </Button>
                    )}
                    {copyState === "failed" && (
                        <span className="opacity-70">Could not copy; logged to the console.</span>
                    )}
                </div>
            )}
            {failedToLoad.length > 0 && (
                <Section title="Layers that could not load">
                    {failedToLoad.map(({ entryId, status, name }) => (
                        <div key={entryId}>
                            <div className="font-medium">{name}</div>
                            <EntryStatusDetails
                                status={{ errors: status.errors, notices: [] }}
                                onRetry={retryEntry && (() => void retryEntry(entryId))}
                            />
                        </div>
                    ))}
                </Section>
            )}
            {withNotices.length > 0 && (
                <Section title="Layer notices">
                    {withNotices.map(({ entryId, status, name }) => (
                        <div key={entryId}>
                            <div className="font-medium">{name}</div>
                            <EntryStatusDetails status={{ errors: [], notices: status.notices }} />
                        </div>
                    ))}
                </Section>
            )}
            <Section title="Renderer">
                <p>{chart.rendererLabel}</p>
                {chart.gpuDescription && <p className="opacity-70">{chart.gpuDescription}</p>}
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
    );
});

export default SpatialStatusPanel;
