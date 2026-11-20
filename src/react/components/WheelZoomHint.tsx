import { useEffect, useRef, useState } from "react";
import { WHEEL_ZOOM_PASSED_EVENT } from "@/lib/wheelZoomGate";

const SHOWN_MS = 1200;

const isMac = typeof navigator !== "undefined" && /Mac|iPhone|iPad/.test(navigator.userAgent);
const MESSAGE = isMac ? "Use ⌘ + scroll to zoom" : "Use Ctrl + scroll to zoom";

/**
 * Says how to zoom when a plain wheel over the chart scrolls the page instead (see
 * `wheelZoomGate`). Place it beside the chart's canvas, inside the chart's positioned content:
 * it listens on its parent element for the event the canvas sends up, and covers the content.
 */
export default function WheelZoomHint() {
    const ref = useRef<HTMLDivElement>(null);
    const [shown, setShown] = useState(false);
    useEffect(() => {
        const host = ref.current?.parentElement;
        if (!host) return;
        let timer: ReturnType<typeof setTimeout> | undefined;
        const show = () => {
            setShown(true);
            clearTimeout(timer);
            timer = setTimeout(() => setShown(false), SHOWN_MS);
        };
        host.addEventListener(WHEEL_ZOOM_PASSED_EVENT, show);
        return () => {
            host.removeEventListener(WHEEL_ZOOM_PASSED_EVENT, show);
            clearTimeout(timer);
        };
    }, []);
    return (
        <div
            ref={ref}
            aria-hidden
            className="pointer-events-none absolute inset-0 z-20 flex items-center justify-center transition-opacity duration-300"
            style={{ opacity: shown ? 1 : 0, backgroundColor: "rgba(15, 23, 42, 0.25)" }}
        >
            <div
                className="rounded-md px-4 py-2 text-sm"
                style={{ backgroundColor: "rgba(15, 23, 42, 0.8)", color: "white" }}
            >
                {MESSAGE}
            </div>
        </div>
    );
}
