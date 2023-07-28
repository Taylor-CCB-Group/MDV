import { io } from "socket.io-client";
import { ChartManager } from "../charts/charts";
import { DataModel } from "../table/DataModel";

declare global {
    interface Window {
        vuplex?: {
            postMessage: (msg: any) => void;
            addEventListener: (type: string, func: (msg: any) => void) => void;
        };
    }
}

type PopoutMessage = {
    type: "popout";
    chartID: string;
}
type FilterMessage = {
    type: "filter";
    dataSource: string;
    indices: Int32Array;
}
type ErrorMessage = {
    type: "error";
    message: string;
}
type MDVMessage = PopoutMessage | FilterMessage | ErrorMessage;

export default async function connectWebsocket(url: string, cm: ChartManager) {
    setupVuplex();
    const socket = io(url);

    function sendMessage(msg: MDVMessage) {
        socket.emit("message", msg);
        if (window.vuplex) window.vuplex.postMessage(msg);
    }

    const originalPopOutChart = cm._popOutChart;
    cm._popOutChart = (chart) => {
        originalPopOutChart.apply(cm, [chart]);
        sendMessage({ type: "popout", chartID: chart.config.id });
    }

    socket.connect();
    function popout(chartID: string) {
        const chart = cm.charts[chartID];
        if (!chart) {
            console.error(`Chart ${chartID} not found`);
        }
        cm._popOutChart(chart.chart);
    }
    socket.on("message", async (chartID: string) => {
        // await new Promise((resolve) => setTimeout(resolve, 5000));
        const chart = cm.charts[chartID];
        if (!chart) {
            console.error(`Chart ${chartID} not found`);
            socket.emit("popout_fail", chartID);
            sendMessage({ type: "error", message: `Chart ${chartID} not found` });
        } else {
            cm._popOutChart(chart.chart);
        }
    });
    // DataModel: register with all datastore changes & send appropriate updates.
    // TOOD: Types, Zod?
    for (const ds of cm.dataSources) {
        const dataModel = new DataModel(ds.dataStore);
        dataModel.addListener('IPC', () => {
            sendMessage({ type: "filter", dataSource: ds.name, indices: dataModel.data });
        });
    }
    
    function setupVuplex() {
        console.log("setupVuplex...");
        function addMessageListener() {
            window.vuplex.postMessage({ type: "vuplex_ready" });
            window.vuplex.addEventListener("message", (msg: string) => {
                const data = JSON.parse(msg) as PopoutMessage;
                console.log("vuplext message", JSON.stringify(msg, null, 2));
                if (data.type === "popout") {
                    console.log("chartID: ", data.chartID);
                    popout(data.chartID);
                }
            });
        }
        if (window.vuplex) {
            console.log("vuplex already exists, adding listener");
            addMessageListener();
        } else {
            console.log("vuplex doesn't exist yet, waiting for vuplexready");
            window.addEventListener('vuplexready', () => {
                console.log("vuplexready event received");
                addMessageListener();
            });
        }
    }
    
    // debugger;
}
