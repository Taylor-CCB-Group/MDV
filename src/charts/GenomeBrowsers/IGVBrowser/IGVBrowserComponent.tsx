import { observer } from "mobx-react-lite";
import { useChart } from "@/react/context";
import { useEffect, useRef, useState } from "react";
import type IGVBrowser from "./IGVBrowser";
import { BASE_TRACK_ID } from "./IGVBrowser";
// @ts-expect-error igv typings are incomplete in this setup
import igv from "igv/dist/igv.esm.js";
import type { Browser } from "igv";
import "./mdvFeatureTrack";



const IGVBrowserComponent = observer(() => {
    const chart = useChart() as IGVBrowser;
    const rootRef = useRef<HTMLDivElement | null>(null);
    const browserRef = useRef<Browser | null>(null);
    const lastLocusRef = useRef<string | null>(null);
    const [loading, setLoading] = useState(true);
    const [searching, setSearching] = useState(false);
    const [error, setError] = useState<string | null>(null);


    useEffect(() => {
        let disposed = false;
        let localBrowser: Browser | null = null;
        //this  is probably more complicated then it needs to be but
        //prevents errors on double mounting
        async function init() {
            if (!rootRef.current) return;
            setLoading(true);
            setError(null);
            try {
                //this will get all the features from the datasource and add them to track config
                //this will only work with small datasources
                const config = await chart.getInitialBrowserConfig();
                if (disposed) return;
                //create the browser instance in the chart div
                const browser = await Promise.race([
                    igv.createBrowser(rootRef.current, config),
                    new Promise<never>((_, reject) => {
                        setTimeout(() => reject(new Error("Timed out while initializing IGV browser.")), 20000);
                    }),
                ]);
                if (disposed) {
                    igv.removeBrowser(browser);
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
                    const frame = Array.isArray(frames) ? frames[0] : frames;
                    if (!frame?.chr) return;
                    const locus = `${frame.chr}:${Math.floor(frame.start)}-${Math.floor(frame.end)}`;
                    lastLocusRef.current = locus;
                    chart.setLocationFromLocus(frame.chr, frame.start, frame.end);
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
                igv.removeBrowser(b);
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
