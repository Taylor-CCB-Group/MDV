/**
 * PROTOTYPE (layout lab) — throwaway. Moves a legacy chart's existing DOM into a card an engine
 * rendered, resizes the chart when the card's size settles, and puts it back on unmount.
 *
 * Hosting needs three things the pane workspace gets in the way of:
 * - its drag (by the title) and resize handles stay on the chart div, so they are removed;
 * - chart overlays have their own stacking, so the slot is isolated and the card's handles win;
 * - a redraw per frame of a gesture is wasted, so it waits until the gesture ends.
 */
import { useLayoutEffect, useRef } from "react";
import type BaseChart from "@/charts/BaseChart";
import { removeDraggable, removeResizable } from "@/utilities/Elements";
import { countSetSize, isGestureActive, onGestureEnd } from "./labShared";

const STYLE_KEYS = ["position", "left", "top", "width", "height"] as const;

export default function ChartHost({ chart }: { chart: BaseChart<any> }) {
    const slotRef = useRef<HTMLDivElement>(null);

    useLayoutEffect(() => {
        const slot = slotRef.current;
        if (!slot) return;
        const div = chart.getDiv();
        const home = div.parentElement;
        const next = div.nextSibling;
        const saved = STYLE_KEYS.map((k) => div.style[k]);

        // the pane workspace's own drag and resize; not put back, as the lab owns the page from here
        removeDraggable(div);
        removeResizable(div);

        div.style.position = "absolute";
        div.style.left = "0";
        div.style.top = "0";
        div.style.width = "100%";
        div.style.height = "100%";
        slot.appendChild(div);

        let timer = 0;
        let pending = false;
        const redraw = () => {
            pending = false;
            countSetSize();
            chart.setSize();
        };
        const ro = new ResizeObserver(() => {
            window.clearTimeout(timer);
            if (isGestureActive()) {
                pending = true;
                return;
            }
            timer = window.setTimeout(redraw, 30);
        });
        ro.observe(slot);
        const offEnd = onGestureEnd(() => {
            if (!pending) return;
            window.clearTimeout(timer);
            // after the card's settle transition
            timer = window.setTimeout(redraw, 220);
        });

        return () => {
            ro.disconnect();
            offEnd();
            window.clearTimeout(timer);
            STYLE_KEYS.forEach((k, i) => {
                div.style[k] = saved[i];
            });
            if (home?.isConnected) home.insertBefore(div, next && next.parentNode === home ? next : null);
        };
    }, [chart]);

    return <div ref={slotRef} className="lab-slot" />;
}
