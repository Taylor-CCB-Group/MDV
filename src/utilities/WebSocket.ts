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


export default async function connectWebsocket(url: string, cm: ChartManager) {
    setupVuplex();
    const socket = io(url);
    const originalPopOutChart = cm._popOutChart;
    cm._popOutChart = (chart) => {
        originalPopOutChart.apply(cm, [chart]);
        socket.emit("popout", chart.config.id);
        if (window.vuplex) window.vuplex.postMessage({ type: "popout", chartID: chart.config.id });
    }

    socket.connect();
    function popout(chartID: string) {
        const chart = cm.charts[chartID];
        if (!chart) {
            console.error(`Chart ${chartID} not found`);
        }
        cm._popOutChart(chart.chart);
    }
    socket.on("popout", async (chartID: string) => {
        // await new Promise((resolve) => setTimeout(resolve, 5000));
        const chart = cm.charts[chartID];
        if (!chart) {
            console.error(`Chart ${chartID} not found`);
            socket.emit("popout_fail", chartID);
        } else {
            cm._popOutChart(chart.chart);
        }
    });
    // DataModel: register with all datastore changes & send appropriate updates.
    // TOOD: Types, Zod?
    for (const ds of cm.dataSources) {
        const dataModel = new DataModel(ds.dataStore);
        dataModel.addListener('socket', () => {
            socket.emit('filter', {dataSource: ds.name, indices: dataModel.data});
        });
    }
    
    function setupVuplex() {
        if (window.vuplex) {
            window.vuplex.addEventListener("message", (msg) => {
                if (msg.type === "popout") {
                    popout(msg.chartID);
                }
            });
        }
    }
    
    // debugger;
}
