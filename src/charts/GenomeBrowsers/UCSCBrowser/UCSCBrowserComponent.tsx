import { observer } from "mobx-react-lite";
import { useChart } from "@/react/context";
import { useEffect, useState, useMemo } from "react";
import { useDebouncedChartSize } from "@/react/hooks";
import type UCSCBrowser from "./UCSCBrowser";
import type { UCSCBrowserConfig} from "./UCSCBrowser";
import { useFieldSpecs } from "@/react/hooks";
import { useHighlightedIndex } from "@/react/selectionHooks";
import {locationFromFieldValues } from "../genomicLocationUtils";
import { useGenomicInfo } from "../genomicLocationUtils";
import {type GenomeLocation,applyViewMargins} from "../genomicLocationUtils";
import { runInAction } from "mobx";


const UCSCBrowserComponent = observer(() => {
    const chart = useChart() as UCSCBrowser;
    const size = useDebouncedChartSize();
    const config = chart.config as UCSCBrowserConfig;
    const [loading, setLoading] = useState(true);
    const [loadError, setLoadError] = useState<string | null>(null);
    const [coordInput, setCoordInput] = useState("");
    const highlightedIndex = useHighlightedIndex();
 


    const {genomicInfo,areColumnsLoaded} = useGenomicInfo();
    const url_proxy = genomicInfo.ucsc_proxy_url || "ucsc_proxy";



    function getHighlightedRegion():GenomeLocation | null {
        if (!areColumnsLoaded || highlightedIndex == null || highlightedIndex === -1) return null;
        const location = locationFromFieldValues(chart.dataStore,highlightedIndex, genomicInfo);
        if (!location) return null;
        return location[0]; //can only cope with one locus
    }

    const handleZoom = (factor: number) => {
        if (!config.location) return;
            const loc= config.location;
            if (!loc) return;
            const currentLength = loc.end - loc.start;
            const center = (loc.start + loc.end) / 2;
            const newLength = Math.round(currentLength * factor);
            const start = Math.max(0, Math.round(center - newLength / 2));
            const end = Math.round(center + newLength / 2);
            config.location = { chr: loc.chr, start, end };
      
    };

    const handlePan = (direction: "left" | "right") => {
        if (!config.location) return;

        const loc = config.location;
        if (!loc) return;
        const currentLength = loc.end - loc.start;
        const shift = Math.round(currentLength * 0.25); // Pan by 25% of current view
        if (direction === "left") {
            const newStart = Math.max(0, loc.start - shift);
            const newEnd = newStart + currentLength;
            config.location = {chr: loc.chr, start: newStart, end: newEnd};
        } else {
            config.location = {chr: loc.chr, start: loc.start + shift, end: loc.end + shift};
        }

    };

    const handleCoordinateSubmit = () => {
        // Parse formats like: chr1:1000-2000 or chr1:1,000-2,000
        const match = coordInput.match(/^(\w+):([0-9,]+)-([0-9,]+)$/);
        if (match) {
            config.location = {
                chr: match[1],
                start: parseInt(match[2].replace(/,/g, "")),
                end: parseInt(match[3].replace(/,/g, ""))
            };
            setCoordInput("");
        }
    };

    // Format current location for display
    const currentLocationText = config.location 
        ? `${config.location.chr}:${config.location.start.toLocaleString()}-${config.location.end.toLocaleString()}`
        : "";   
 
    const trimmedSrc = config.src?.trim();
    const missingSrcMessage = !trimmedSrc
        ? "Enter a UCSC session URL to load the browser image."
        : "Enter a valid UCSC session URL to load the browser image.";
    
    
    // Calculate src reactively - this will re-run when config.location changes
    const src = useMemo(() => {
        const baseUrl = config.src?.trim();
        if (!baseUrl) {
            return null;
        }

        let url: URL;
        try {
            url = new URL(baseUrl);
        } catch {
            return null;
        }

        // Update or add the position parameter
        if (config.location) {
            const position = `${config.location.chr}:${config.location.start}-${config.location.end}`;
            url.searchParams.set("position", position);
            url.searchParams.set("pix", size[0].toString());
            //highlight the actual region selected
            if (config.highlight_selected_region && highlightedIndex != null){ 
                const region = getHighlightedRegion();
                if (region){
                     const db = url.searchParams.get("db") || "hg38";
                    url.searchParams.set("highlight", `${db}.${region.chr}:${region.start}-${region.end}`);
                }
            }
        }
        // Keep render requests on the same UCSC host as the session URL.
        url.searchParams.set("ucscHost", url.host);

        // Build the proxy URL with all params
        return `${url_proxy}?${url.searchParams.toString()}`;
    }, [trimmedSrc, config.location, size, config.highlight_selected_region, highlightedIndex]);

    // Reset loading/error state when src changes
    useEffect(() => {
        setLoadError(null);
        setLoading(Boolean(src));
    }, [src]);
    //update the location of the browser
    useEffect(() =>{
        if (highlightedIndex !== -1){
            const region = getHighlightedRegion();  
            if (!region) return;
            //run in mobx action to avoid warnings about updating observable state outside of an action
            runInAction(() => {
                config.location = applyViewMargins(region, config.view_margins!);    
            });
        }
    }, [config.view_margins, highlightedIndex]);

    return (
        <div style={{ position: "relative", width: "100%", height: "100%", display: "flex", flexDirection: "column" }}>
            {/* Navigation controls */}
            <div style={{ 
                display: "flex", 
                gap: "4px", 
                padding: "4px 8px", 
                borderBottom: "1px solid #ddd",
                alignItems: "center",
                flexShrink: 0
            }}>
                <button 
                    onClick={() => handlePan("left")} 
                    style={{ padding: "2px 8px", fontSize: "12px" }}
                    title="Pan left"
                >
                    ◀
                </button>
                <button 
                    onClick={() => handleZoom(1.5)} 
                    style={{ padding: "2px 8px", fontSize: "12px" }}
                    title="Zoom out"
                >
                    −
                </button>
                <button 
                    onClick={() => handleZoom(0.67)} 
                    style={{ padding: "2px 8px", fontSize: "12px" }}
                    title="Zoom in"
                >
                    +
                </button>
                <button 
                    onClick={() => handlePan("right")} 
                    style={{ padding: "2px 8px", fontSize: "12px" }}
                    title="Pan right"
                >
                    ▶
                </button>
                <input
                    type="text"
                    value={coordInput || currentLocationText}
                    onChange={(e) => setCoordInput(e.target.value)}
                    onKeyDown={(e) => e.key === "Enter" && handleCoordinateSubmit()}
                    placeholder="chr1:1000-2000"
                    style={{ 
                        marginLeft: "8px", 
                        padding: "2px 6px", 
                        fontSize: "12px",
                        width: "150px"
                    }}
                />
                <button 
                    onClick={handleCoordinateSubmit}
                    style={{ padding: "2px 8px", fontSize: "12px" }}
                >
                    Go
                </button>
            </div>
            
            {/* Image container */}
            <div style={{ position: "relative", flex: 1, overflow: "auto" }}>
                {loading && !loadError && (
                    <div style={{
                        position: "absolute",
                        top: 0,
                        left: 0,
                        right: 0,
                        bottom: 0,
                        display: "flex",
                        alignItems: "center",
                        justifyContent: "center",
                        backgroundColor: "rgba(255, 255, 255, 0.8)",
                        zIndex: 1
                    }}>
                        <div>Loading...</div>
                    </div>
                )}
                {!src ? (
                    <div style={{
                        display: "flex",
                        alignItems: "center",
                        justifyContent: "center",
                        height: "100%",
                        padding: "16px",
                        color: "#666",
                        textAlign: "center"
                    }}>
                        {missingSrcMessage}
                    </div>
                ) : loadError ? (
                    <div style={{
                        display: "flex",
                        alignItems: "center",
                        justifyContent: "center",
                        height: "100%",
                        padding: "16px",
                        color: "#b42318",
                        textAlign: "center"
                    }}>
                        {loadError}
                    </div>
                ) : (
                    <img
                        src={src}
                        title="UCSC Genome Browser"
                        onLoad={() => setLoading(false)}
                        onError={() => {
                            setLoading(false);
                            setLoadError("Unable to load the UCSC browser image.");
                        }}
                    />
                )}
            </div>
        </div>
    );
});


export default UCSCBrowserComponent;
