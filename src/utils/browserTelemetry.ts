/**
 * Browser telemetry with Grafana Faro, off unless the server switches it on.
 *
 * The page asks `telemetry/config` on its own deployment. When the server has
 * `FARO_COLLECTOR_URL` set, the page loads Faro, which records page loads,
 * Web Vitals, how long each request takes and uncaught errors, and posts them to
 * `telemetry/faro` on the same deployment. The server forwards them to Alloy
 * (ADR 0023). Faro's tracing also adds a `traceparent` header to MDV's own
 * requests, so a slow request in the browser links to its trace on the server.
 *
 * Faro is only downloaded when telemetry is on. A page whose API is on another
 * origin, such as a hosted preview pointed at a server with `?dir=`, never
 * starts it, because the server only accepts posts from its own pages.
 */
import { buildApiUrl } from "@/utils/mdvRouting";

type TelemetryConfig = { enabled: false } | { enabled: true; app_name: string; sample_rate: number };

let started = false;

export async function startBrowserTelemetry() {
    if (started) return;
    started = true;
    try {
        const collectUrl = new URL(buildApiUrl("telemetry/faro"), window.location.href);
        if (collectUrl.origin !== window.location.origin) return;

        const response = await fetch(buildApiUrl("telemetry/config"));
        if (!response.ok) return;
        const config: TelemetryConfig = await response.json();
        if (!config.enabled) return;

        const [{ initializeFaro, getWebInstrumentations }, { TracingInstrumentation }] = await Promise.all([
            import("@grafana/faro-web-sdk"),
            import("@grafana/faro-web-tracing"),
        ]);
        initializeFaro({
            url: collectUrl.toString(),
            app: {
                name: config.app_name,
                version: import.meta.env.VITE_GIT_COMMIT_HASH || "dev",
            },
            sessionTracking: { samplingRate: config.sample_rate },
            // MDV writes a lot to the console. Uncaught errors are still recorded.
            instrumentations: [...getWebInstrumentations({ captureConsole: false }), new TracingInstrumentation()],
        });
    } catch (error) {
        // Telemetry must never get in the way of using MDV.
        console.warn("Browser telemetry did not start", error);
    }
}
