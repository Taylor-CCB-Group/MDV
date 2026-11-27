/**
 * PROTOTYPE (layout lab) — throwaway, dev builds only. Open a project with `?layoutLab=inhouse` or
 * `?layoutLab=gridstack` to lay its charts out in one view with the in-house engine or gridstack 14.
 * Not for main: the answer gets recorded in an ADR and this folder is deleted.
 */
import { createElement } from "react";
import type ChartManager from "@/charts/ChartManager";
import { createMdvPortal } from "@/react/react_utils";
import LayoutLab, { type LabMode } from "./LayoutLab";

export function startLayoutLab(cm: ChartManager, mode: LabMode) {
    const mount = () => {
        // a sibling over the content area: the content area itself is cleared on view changes
        const host = document.createElement("div");
        cm.containerDiv.style.position = "relative";
        Object.assign(host.style, {
            position: "absolute",
            left: "0",
            right: "0",
            bottom: "0",
            top: `${cm.contentDiv.offsetTop}px`,
        });
        cm.containerDiv.appendChild(host);
        createMdvPortal(createElement(LayoutLab, { cm, initialMode: mode }), host, undefined, cm);
    };
    // wait for the first view so all its charts exist (`_toLoadCharts` is cleared once they have)
    const loaded = (cm as unknown as { _toLoadCharts?: Set<unknown> })._toLoadCharts === undefined;
    if (loaded && Object.keys(cm.charts).length) {
        mount();
        return;
    }
    const id = "layout-lab-start";
    cm.addListener(id, (type: string) => {
        if (type !== "view_loaded") return;
        cm.removeListener(id);
        mount();
    });
}
