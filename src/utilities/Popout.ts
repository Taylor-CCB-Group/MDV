import type { Chart } from "@/charts/charts";
import { removeDraggable, removeResizable } from "./Elements";

export default function popoutChart(chart: Chart) {
    const { chartManager } = window.mdv;
    const mainWindow = window;
    const div = chart.getDiv();
    const originalParent = div.parentElement;
    const wtop = window.screenTop ? window.screenTop : window.screenY;
    const wleft = window.screenLeft ? window.screenLeft : window.screenX;
    const { width, height, left, top } = div.getBoundingClientRect();
    const w = Math.floor(Math.max(100, width));
    const h = Math.floor(Math.max(100, height));
    const l = Math.floor(left + wleft);
    const t = Math.floor(top + wtop);
    if (div.gridstackPopoutCallback) div.gridstackPopoutCallback();
    removeResizable(div);
    removeDraggable(div);
    const popoutWindow = window.open(
        "",
        chart.config.title,
        // "width=800,height=600",
        `width=${w},height=${h},left=${l},top=${t}`,
    );
    popoutWindow.document.body.style.overflow = "hidden";

    // Function to add a stylesheet or style element to the popout window
    const addStyleElement = (
        styleElement: HTMLStyleElement | HTMLLinkElement,
    ) => {
        if (styleElement instanceof HTMLLinkElement) {
            const newLink = popoutWindow.document.createElement("link");
            newLink.rel = "stylesheet";
            newLink.href = styleElement.href;
            popoutWindow.document.head.appendChild(newLink);
        } else if (styleElement instanceof HTMLStyleElement) {
            const newStyle = popoutWindow.document.createElement("style");
            newStyle.innerHTML = styleElement.innerHTML;
            popoutWindow.document.head.appendChild(newStyle);
        }
    };

    // Initial copy of all existing styles
    Array.from(
        document.querySelectorAll('link[rel="stylesheet"], style'),
    ).forEach((styleElement) => {
        addStyleElement(styleElement as HTMLStyleElement | HTMLLinkElement);
    });

    // Observe the main window's head for new stylesheets and style elements
    const observer = new MutationObserver((mutations) => {
        mutations.forEach((mutation) => {
            if (
                mutation.type === "childList" &&
                mutation.addedNodes.length > 0
            ) {
                mutation.addedNodes.forEach((node) => {
                    if (
                        node instanceof HTMLLinkElement ||
                        node instanceof HTMLStyleElement
                    ) {
                        addStyleElement(node);
                    }
                });
            }
        });
    });

    // Start observing the main window's head for changes
    observer.observe(document.head, {
        childList: true,
        subtree: true,
    });

    popoutWindow.document.body.append(div);
    chart.changeBaseDocument(popoutWindow.document);
    const popStyles = pushSetStyles(chart);

    const body = popoutWindow.document.body;
    popoutWindow.addEventListener("resize", () => {
        chart.setSize(body.clientWidth, body.clientHeight);
    });

    // I don't think we use this win property much, but it's here for consistency
    // with previous implementation & comments
    chartManager.charts[chart.config.id].win = popoutWindow;
    popoutWindow.addEventListener("beforeunload", () => {
        popStyles();
        originalParent.appendChild(div);
        chartManager._makeChartRD(chart);
        chartManager.charts[chart.config.id].win = mainWindow;
        chart.changeBaseDocument(document);
        observer.disconnect();
    });
    mainWindow.addEventListener("unload", () => {
        popoutWindow.close();
    });
    return popoutWindow;
}

/** apply the styles for popped-out version of chart & return a function that will restore them later */
function pushSetStyles(chart: Chart) {
    const div = chart.getDiv();
    const styles = {
        height: div.style.height,
        width: div.style.width,
        top: div.style.top,
        left: div.style.left,
        iconDisplay: chart.popoutIcon.style.display,
    };
    div.style.height = "100vh";
    div.style.width = "100vw";
    div.style.top = "0";
    div.style.left = "0";
    chart.popoutIcon.style.display = "none";
    return () => {
        chart.popoutIcon.style.display = styles.iconDisplay;
        div.style.height = styles.height;
        div.style.width = styles.width;
        div.style.top = styles.top;
        div.style.left = styles.left;
        chart.setSize();
    };
}
