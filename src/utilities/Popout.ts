import { Chart } from "@/charts/charts";
import { removeDraggable, removeResizable } from "./Elements";

export default function popoutChart(chart: Chart) {
    const { chartManager } = window.mdv;
    const mainWindow = window;
    const div = chart.getDiv();
    const originalParent = div.parentElement;
    if (div.gridstackPopoutCallback) div.gridstackPopoutCallback();
    removeResizable(div);
    removeDraggable(div);
    const popoutWindow = window.open('', chart.config.title, 'width=800,height=600');
    popoutWindow.document.body.style.overflow = 'hidden';
    Array.from(document.styleSheets).forEach(styleSheet => {
        if (styleSheet.href) {
            const newLink = popoutWindow.document.createElement('link');
            newLink.rel = 'stylesheet';
            newLink.href = styleSheet.href;
            popoutWindow.document.head.appendChild(newLink);
        } else if (styleSheet.cssRules) {
            const newStyle = popoutWindow.document.createElement('style');
            Array.from(styleSheet.cssRules).forEach(rule => {
                newStyle.appendChild(popoutWindow.document.createTextNode(rule.cssText));
            });
            popoutWindow.document.head.appendChild(newStyle);
        }
    });
    popoutWindow.document.body.append(div);
    chart.changeBaseDocument(popoutWindow.document);
    const popStyles = pushSetStyles(chart);

    const body = popoutWindow.document.body;
    const resizeObserver = new ResizeObserver(() => {
        chart.setSize(body.clientWidth, body.clientHeight);
    }); 
    resizeObserver.observe(popoutWindow.document.documentElement);

    // I don't think we use this win property much, but it's here for consistency
    // with previous implementation & comments
    chartManager.charts[chart.config.id].win = popoutWindow;
    popoutWindow.addEventListener('beforeunload', () => {
        popStyles();
        originalParent.appendChild(div);
        chartManager._makeChartRD(chart);
        chartManager.charts[chart.config.id].win = mainWindow;
        chart.changeBaseDocument(document);
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
    div.style.height = '100vh';
    div.style.top = '0';
    div.style.left = '0';
    chart.popoutIcon.style.display = 'none';
    return () => {
        chart.popoutIcon.style.display = styles.iconDisplay;
        div.style.height = styles.height;
        div.style.width = styles.width;
        div.style.top = styles.top;
        div.style.left = styles.left;
        chart.setSize();
    }
}
