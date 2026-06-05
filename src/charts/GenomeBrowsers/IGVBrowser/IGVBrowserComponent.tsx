import { observer } from "mobx-react-lite";
import { useChart } from "@/react/context";
import { useEffect, useRef, useState } from "react";
import type IGVBrowser from "./IGVBrowser";
import { BASE_TRACK_ID } from "./IGVBrowser";
import type { Browser } from "igv";
import { useFieldSpec } from "@/react/hooks";




const IGV_DOCUMENT_EVENTS = new Set([
    "mousemove",
    "mouseup",
    "mouseleave",
    "mouseexit",
    "keyup",
]);

const IGV_WINDOW_EVENTS = new Set([
    "resize",
    "beforeunload",
]);

let igvEventRerouteCleanup: (() => void) | null = null;
let igvEventRerouteRefCount = 0;

function installIgvGlobalEventReroute(targetDoc: Document) {
    const targetWin = targetDoc.defaultView;
    if (!targetWin || targetDoc === document) {
        return () => {};
    }

    if (igvEventRerouteCleanup) {
        igvEventRerouteRefCount += 1;
        return () => {
            igvEventRerouteRefCount -= 1;
            if (igvEventRerouteRefCount <= 0) {
                igvEventRerouteCleanup?.();
                igvEventRerouteCleanup = null;
                igvEventRerouteRefCount = 0;
            }
        };
    }

    igvEventRerouteRefCount = 1;

    const originalDocumentAdd = document.addEventListener.bind(document);
    const originalDocumentRemove = document.removeEventListener.bind(document);
    const originalWindowAdd = window.addEventListener.bind(window);
    const originalWindowRemove = window.removeEventListener.bind(window);

    document.addEventListener = function (
        type: string,
        listener: EventListenerOrEventListenerObject | null,
        options?: boolean | AddEventListenerOptions,
    ) {
        if (!listener) {
            return;
        }
        if (IGV_DOCUMENT_EVENTS.has(type)) {
            return targetDoc.addEventListener(type, listener, options);
        }
        return originalDocumentAdd(type, listener, options);
    };

    document.removeEventListener = function (
        type: string,
        listener: EventListenerOrEventListenerObject | null,
        options?: boolean | EventListenerOptions,
    ) {
        if (!listener) {
            return;
        }
        if (IGV_DOCUMENT_EVENTS.has(type)) {
            return targetDoc.removeEventListener(type, listener, options);
        }
        return originalDocumentRemove(type, listener, options);
    };

    window.addEventListener = function (
        type: string,
        listener: EventListenerOrEventListenerObject | null,
        options?: boolean | AddEventListenerOptions,
    ) {
        if (!listener) {
            return;
        }
        if (IGV_WINDOW_EVENTS.has(type)) {
            return targetWin.addEventListener(type, listener, options);
        }
        return originalWindowAdd(type, listener, options);
    };

    window.removeEventListener = function (
        type: string,
        listener: EventListenerOrEventListenerObject | null,
        options?: boolean | EventListenerOptions,
    ) {
        if (!listener) {
            return;
        }
        if (IGV_WINDOW_EVENTS.has(type)) {
            return targetWin.removeEventListener(type, listener, options);
        }
        return originalWindowRemove(type, listener, options);
    };

    igvEventRerouteCleanup = () => {
        document.addEventListener = originalDocumentAdd;
        document.removeEventListener = originalDocumentRemove;
        window.addEventListener = originalWindowAdd;
        window.removeEventListener = originalWindowRemove;
    };

    return () => {
        igvEventRerouteRefCount -= 1;
        if (igvEventRerouteRefCount <= 0) {
            igvEventRerouteCleanup?.();
            igvEventRerouteCleanup = null;
            igvEventRerouteRefCount = 0;
        }
    };
}




const IGVBrowserComponent = observer(() => {
    const chart = useChart() as IGVBrowser;
    const rootRef = useRef<HTMLDivElement | null>(null);
    const browserRef = useRef<Browser | null>(null);
    const [loading, setLoading] = useState(true);
    const [searching, setSearching] = useState(false);
    const [error, setError] = useState<string | null>(null);
    const featureLabelColumn = useFieldSpec(chart.config.feature_label==="_none_"?undefined:chart.config.feature_label);
    const remountNonce = chart.browserRemountNonce;

    useEffect(() => {
        return installIgvGlobalEventReroute(chart.__doc__);
    }, [chart, remountNonce]);

    useEffect(() => {
        let disposed = false;
        let localBrowser: Browser | null = null;
        let igvApi: any = null;
        //this  is probably more complicated then it needs to be but
        //prevents errors on double mounting
        async function init() {
            if (!rootRef.current) return;
            setLoading(true);
            setError(null);
            try {
                await import("./mdvFeatureTrack");
                const igvModule = await import("igv");
                igvApi = igvModule.default;
                //this will get all the features from the datasource and add them to track config
                //this will only work with small datasources
                const config = await chart.getInitialBrowserConfig();
                if (disposed) return;
                //create the browser instance in the chart div
                const browser = await Promise.race([
                    igvApi.createBrowser(rootRef.current, config),
                    new Promise<never>((_, reject) => {
                        setTimeout(() => reject(new Error("Timed out while initializing IGV browser.")), 30000);
                    }),
                ]);
                if (disposed) {
                    igvApi?.removeBrowser(browser);
                    return;
                }
                //now the browser has loaded - keep the its shadow dom for popout
                chart.keepShadowDom(rootRef.current);
                localBrowser = browser;
                browserRef.current = browser;
                const baseTrack = browser.findTracks("id", BASE_TRACK_ID)?.[0];
                chart.attachBrowser(browser, baseTrack);

                browser.on("trackclick", (track: any, popupData: any) => {
                    if (track?.id !== BASE_TRACK_ID && track?.config?.id !== BASE_TRACK_ID) {
                        return undefined;
                    }
                    const id = popupData?.[0]?.id;
                    if (id ){
                        chart.dataStore.dataHighlighted([Number.parseInt(id)], chart);
                    }
                    return false;
                });

                browser.on("locuschange", (frames: any) => {
                    //const frame = Array.isArray(frames) ? frames[0] : frames;
                    //if (!frame?.chr) return;
                    //const locus = `${frame.chr}:${Math.floor(frame.start)}-${Math.floor(frame.end)}`;
                    //lastLocusRef.current = locus;
                    //chart.setLocationFromLocus(frame.chr, frame.start, frame.end);
                });
                setLoading(false);
            } catch (e) {
                if (!disposed) {
                    setError((e as Error).message || "Failed to initialize IGV browser");
                }
            } finally {
                if (!disposed) {
                    setLoading(false);
                }
            }
        }

        init();

        return () => {
            disposed = true;
            const b = localBrowser || browserRef.current;
            if (b) {
                igvApi?.removeBrowser(b);
            }
            if (browserRef.current === b) browserRef.current = null;
        };
    }, []);

    useEffect(() => {
        chart.setSearchPendingHandler(setSearching);
        return () => {
            chart.setSearchPendingHandler(null);
        };
    }, [chart]);

    useEffect(() => {
        const lf = chart.config.feature_label;
        if (lf !== "_none_" &&(!featureLabelColumn || featureLabelColumn.field !== lf)) return;
        void chart.updateMDVFeatures();
    }, [chart, chart.config.feature_label, featureLabelColumn, ]);

     /*useEffect(() => {
        const locus = buildLocusSearch(chart.config.location);
        if (!locus) return;
        if (locus === lastLocusRef.current) return;
        browserRef.current?.search(locus);
    }, [chart.config.location]);*/

    /*useEffect(() => {
        if (highlightedIndex >= 0) {
            chart.onDataHighlighted({ indexes: [highlightedIndex] });
        }
    }, [highlightedIndex, chart]);*/

    return (
        <div style={{ width: "100%", height: "100%", position: "relative" }}>
            {loading && !browserRef.current ? (
                <div style={{ position: "absolute", inset: 0, display: "flex", alignItems: "center", justifyContent: "center", zIndex: 1, background: "rgba(255,255,255,0.6)", pointerEvents: "none" }}>
                    Loading IGV...
                </div>
            ) : null}
            {searching && browserRef.current ? (
                <div style={{ position: "absolute", inset: 0, display: "flex", alignItems: "center", justifyContent: "center", zIndex: 2, background: "rgba(255,255,255,0.35)", pointerEvents: "none" }}>
                    <div style={{ display: "flex", alignItems: "center", gap: "8px", background: "rgba(255,255,255,0.9)", border: "1px solid #d1d5db", borderRadius: "9999px", padding: "6px 12px", fontSize: "12px", color: "#1f2937" }}>
                        <span
                            style={{
                                width: "12px",
                                height: "12px",
                                border: "2px solid #cbd5e1",
                                borderTopColor: "#2563eb",
                                borderRadius: "9999px",
                                display: "inline-block",
                                animation: "mdv-igv-spin 0.8s linear infinite",
                            }}
                        />
                        Searching...
                    </div>
                </div>
            ) : null}
            <style>{`@keyframes mdv-igv-spin { from { transform: rotate(0deg); } to { transform: rotate(360deg); } }`}</style>
            {error ? (
                <div style={{ color: "#b91c1c", padding: "8px", fontSize: "12px" }}>{error}</div>
            ) : null}
            {chart.guardrailMessage ? (
                <div style={{ color: "#92400e", background: "#fffbeb", border: "1px solid #fde68a", borderRadius: "4px", margin: "4px", padding: "4px 8px", fontSize: "11px" }}>
                    {chart.guardrailMessage}
                </div>
            ) : null}
            <div ref={rootRef} style={{ width: "100%", height: "100%", overflow: "auto" }} />
        </div>
    );
});

export default IGVBrowserComponent;
