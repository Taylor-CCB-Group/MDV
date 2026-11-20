import { useEffect, useRef, useState } from "react";
import {
    hasLearnedViewDrag,
    VIEW_DRAG_PASSED_EVENT,
    VIEW_MODIFIER_NAME,
    WHEEL_ZOOM_PASSED_EVENT,
} from "@/lib/viewGestureGate";

const SHOWN_MS = 1200;

const WHEEL_MESSAGE = `Use ${VIEW_MODIFIER_NAME} + scroll to zoom`;
const DRAG_MESSAGE = `Hold ${VIEW_MODIFIER_NAME} and drag to move the view`;

/**
 * Says how to zoom or move the view when a plain gesture over the chart doesn't (see
 * `viewGestureGate`): on every plain wheel that scrolls the page instead, and on plain drags that
 * draw a selection instead, until a drag has moved a view in this browser. Place it beside the
 * chart's canvas, inside the chart's positioned content: it listens on its parent element for the
 * events the canvas sends up, and covers the content.
 */
export default function ViewGestureHint() {
    const ref = useRef<HTMLDivElement>(null);
    const [message, setMessage] = useState<string | null>(null);
    useEffect(() => {
        const host = ref.current?.parentElement;
        if (!host) return;
        let timer: ReturnType<typeof setTimeout> | undefined;
        const show = (text: string) => {
            setMessage(text);
            clearTimeout(timer);
            timer = setTimeout(() => setMessage(null), SHOWN_MS);
        };
        const onWheel = () => show(WHEEL_MESSAGE);
        const onDrag = () => {
            if (!hasLearnedViewDrag()) show(DRAG_MESSAGE);
        };
        host.addEventListener(WHEEL_ZOOM_PASSED_EVENT, onWheel);
        host.addEventListener(VIEW_DRAG_PASSED_EVENT, onDrag);
        return () => {
            host.removeEventListener(WHEEL_ZOOM_PASSED_EVENT, onWheel);
            host.removeEventListener(VIEW_DRAG_PASSED_EVENT, onDrag);
            clearTimeout(timer);
        };
    }, []);
    return (
        <div
            ref={ref}
            aria-hidden
            data-view-gesture-hint
            className="pointer-events-none absolute inset-0 z-20 flex items-center justify-center transition-opacity duration-300"
            style={{ opacity: message ? 1 : 0, backgroundColor: "rgba(15, 23, 42, 0.25)" }}
        >
            <div
                className="rounded-md px-4 py-2 text-sm"
                style={{ backgroundColor: "rgba(15, 23, 42, 0.8)", color: "white" }}
            >
                {message}
            </div>
        </div>
    );
}
