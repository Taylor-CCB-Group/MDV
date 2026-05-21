import { observer } from "mobx-react-lite";
import { useChart } from "@/react/context";
import { useEffect, useState, useMemo } from "react";
import { useDebouncedChartSize } from "@/react/hooks";
import type UCSCBrowser from "./UCSCBrowser";
import { getLocation,UCSCBrowserConfig,UCSCBrowserLocation,UCSCBrowserViewMargins } from "./UCSCBrowser";
import { useFieldSpecs } from "@/react/hooks";
import { useHighlightedIndex } from "@/react/selectionHooks";
import { highlightedIndexToLocation } from "../genomicLocationUtils";


const UCSCBrowserComponent = observer(() => {
    const chart = useChart() as UCSCBrowser;
    const size = useDebouncedChartSize();
    const config = chart.config as UCSCBrowserConfig;
    const [loading, setLoading] = useState(true);
    const [coordInput, setCoordInput] = useState("");
    const fieldSpecs = useFieldSpecs(chart.config.param);
    const highlightedIndex = useHighlightedIndex();
    const genome = chart.dataStore?.genome;
    const url_proxy = chart.dataStore?.genome?.ucsc_proxy_url || "ucsc_proxy";

    function getHighlightedRegion(): UCSCBrowserLocation | null {
        const location = highlightedIndexToLocation(highlightedIndex, fieldSpecs, Boolean(genome?.svs));
        if (!location) return null;
        return Array.isArray(location) ? location[0] : location;
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
 
    
    // Calculate src reactively - this will re-run when config.location changes
    const src = useMemo(() => {
        // Use the base UCSC URL or config.src
        const baseUrl = config.src || `https://genome.ucsc.edu/cgi-bin/hgTracks?db=${genome?.assembly || "hg38"}`;
        const url = new URL(baseUrl);

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
    }, [config.src, config.location, size, config.highlight_selected_region, highlightedIndex, fieldSpecs, genome?.svs, url_proxy]);

    // Set loading to true when src changes
    useEffect(() => {
        setLoading(true);
    }, [src]);
    //update the location of the browser
    useEffect(() =>{
        if (highlightedIndex !== -1){
            const region = getHighlightedRegion();  
            if (!region) return;
            config.location = getLocation(region, config.view_margins!);    
        }
    }, [config.view_margins, highlightedIndex, fieldSpecs, genome?.svs]);

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
                {loading && (
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
                <img 
                    src={src} 
                    title="UCSC Genome Browser"
                    onLoad={() => setLoading(false)}
                    onError={() => setLoading(false)}
                />
            </div>
        </div>
    );
});


export type { UCSCBrowserConfig, UCSCBrowserLocation,UCSCBrowserViewMargins };
export default UCSCBrowserComponent;
