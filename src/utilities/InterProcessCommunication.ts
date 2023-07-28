import { io } from "socket.io-client";
import { ChartManager } from "../charts/charts";
import { DataModel } from "../table/DataModel";

declare global {
    interface Window {
        vuplex?: {
            postMessage: (msg: any) => void;
            addEventListener: (type: string, func: (event: {data: MDVMessage}) => void) => void;
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
            console.log("current chart keys: ", Object.keys(cm.charts).join(", "));
            sendMessage({ type: "error", message: `Popout failed: chart ${chartID} not found... trying again...` });
            // maybe introduce timeout retry here... not sure why not all charts are available at the same time.
            // there should be a different mechanism for this, but for now, this will do.
            setTimeout(() => popout(chartID), 100);
        } else {
            cm._popOutChart(chart.chart);
        }
    }
    // should msg be a string? or a JSON object?
    socket.on("message", async (msg: string) => {
        const data = JSON.parse(msg) as MDVMessage;
        if (data.type === "popout") {
            const { chartID } = data;
            popout(chartID);
        }
    });
    // DataModel: register with all datastore changes & send appropriate updates.
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
            window.vuplex.addEventListener("message", (e) => {
                let msg = e.data;
                if (typeof msg === "string") msg = JSON.parse(msg) as MDVMessage;
                console.log("vuplex message", JSON.stringify(msg, null, 2));
                // TOOD: Types, Zod?
                try {
                    const data = msg;//JSON.parse(msg) as MDVMessage;
                    if (data.type === "popout") {
                        console.log("chartID: ", data.chartID);
                        popout(data.chartID);
                    }
                } catch (e) {
                    console.error("Error parsing vuplex message", e);
                    sendMessage({ type: "error", message: `Error parsing vuplex message: ${e}` });
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
