import { observer } from "mobx-react-lite";
import { useEffect, useRef, useState } from "react";
import { chartDragSetting } from "@/lib/browserSettings";
import {
    DRAG_CHOICE_EVENT,
    hasLearnedViewDrag,
    VIEW_DRAG_PASSED_EVENT,
    VIEW_MODIFIER_NAME,
    WHEEL_ZOOM_PASSED_EVENT,
} from "@/lib/viewGestureGate";

const SHOWN_MS = 1200;

const WHEEL_MESSAGE = `Use ${VIEW_MODIFIER_NAME} + scroll to zoom`;
const DRAG_MESSAGE = `Hold ${VIEW_MODIFIER_NAME} and drag to move the view`;

const PANEL_STYLE = { backgroundColor: "rgba(15, 23, 42, 0.85)", color: "white" };

/**
 * Asks what a plain drag should do, the first time one moves a chart's view in this browser
 * (`chartDragSetting`). It stays until answered; the answer applies to every chart whose tool
 * hasn't been picked, now and from then on.
 */
const DragChoice = observer(function DragChoice({ onDone }: { onDone: () => void }) {
    if (chartDragSetting.isChosen()) return null;
    const choose = (value: "view" | "select") => {
        chartDragSetting.set(value);
        onDone();
    };
    const buttonClass = "rounded-lg px-3 py-1 text-sm font-medium bg-white/15 hover:bg-white/30";
    return (
        <div
            role="dialog"
            aria-label="What dragging does on charts"
            data-drag-choice
            className="absolute bottom-3 left-1/2 z-30 w-max max-w-[90%] -translate-x-1/2 rounded-xl px-4 py-3 text-center text-sm shadow-lg"
            style={PANEL_STYLE}
        >
            <div className="mb-2">Dragging can move the view or draw a selection.</div>
            <div className="mb-2 flex justify-center gap-2">
                <button type="button" className={buttonClass} onClick={() => choose("view")}>
                    Move the view
                </button>
                <button type="button" className={buttonClass} onClick={() => choose("select")}>
                    Draw a selection
                </button>
            </div>
            <div className="text-xs opacity-80">
                Either way, {VIEW_MODIFIER_NAME}-drag moves the view. You can change this later in Settings, under
                Preferences.
            </div>
        </div>
    );
});

/**
 * Says how to zoom or move the view when a plain gesture over the chart doesn't (see
 * `viewGestureGate`): on every plain wheel that scrolls the page instead, and on plain drags that
 * draw a selection instead, until a drag has moved a view in this browser. The first plain drag
 * that moves the view, before this browser has chosen what dragging does, asks (`DragChoice`).
 * Place it beside the chart's canvas, inside the chart's positioned content: it listens on its
 * parent element for the events the canvas sends up, and covers the content.
 */
export default function ViewGestureHint() {
    const ref = useRef<HTMLDivElement>(null);
    const [message, setMessage] = useState<string | null>(null);
    const [asking, setAsking] = useState(false);
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
        const onChoice = () => setAsking(true);
        host.addEventListener(WHEEL_ZOOM_PASSED_EVENT, onWheel);
        host.addEventListener(VIEW_DRAG_PASSED_EVENT, onDrag);
        host.addEventListener(DRAG_CHOICE_EVENT, onChoice);
        return () => {
            host.removeEventListener(WHEEL_ZOOM_PASSED_EVENT, onWheel);
            host.removeEventListener(VIEW_DRAG_PASSED_EVENT, onDrag);
            host.removeEventListener(DRAG_CHOICE_EVENT, onChoice);
            clearTimeout(timer);
        };
    }, []);
    return (
        <>
            <div
                ref={ref}
                aria-hidden
                data-view-gesture-hint
                className="pointer-events-none absolute inset-0 z-20 flex items-center justify-center transition-opacity duration-300"
                style={{ opacity: message ? 1 : 0, backgroundColor: "rgba(15, 23, 42, 0.25)" }}
            >
                <div className="rounded-xl px-4 py-2 text-center text-sm" style={PANEL_STYLE}>
                    {message}
                </div>
            </div>
            {asking && <DragChoice onDone={() => setAsking(false)} />}
        </>
    );
}
