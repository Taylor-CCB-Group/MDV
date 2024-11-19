import { io } from "socket.io-client";
import type ChartManager from "../charts/ChartManager";
import { DataModel } from "../table/DataModel";

type VuplexCallback = (event: { data: MDVMessage }) => void;

declare global {
    interface Window {
        vuplex?: {
            postMessage: (msg: any) => void;
            addEventListener: (type: string, func: VuplexCallback) => void;
        };
    }
}

type PopoutMessage = {
    type: "popout";
    chartID: string;
};
type FilterMessage = {
    type: "filter";
    dataSource: string;
    indices: Int32Array;
};
type ErrorMessage = {
    type: "error";
    message: string;
};
type ChartEventNames =
    | "state_saved"
    | "chart_added"
    | "chart_removed"
    | "view_loaded";
type Chart = unknown; //todo
type ViewId = string;
type ChartState = ReturnType<typeof ChartManager.prototype.getState>;
type ChartMsgValue<T extends ChartEventNames> = T extends
    | "chart_added"
    | "chart_removed"
    ? Chart
    : T extends "view_loaded"
      ? ViewId
      : T extends "state_saved"
        ? { state: ChartState }
        : never;

type MDVMessage = PopoutMessage | FilterMessage | ErrorMessage;

/** initial experimental IPC support. As of this writing, will attempt to
 * - connect a 'vuplex' PostMessage interface for embedding in a VR environment in Unity
 * - connect a websocket to a socket.io server that will forward messages to other clients (not used yet)
 */
export default async function connectIPC(cm: ChartManager) {
    const params = new URLSearchParams(window.location.search);
    const url = params.get("socket") || undefined;
    let initialized = false;
    cm.addListener(
        "ipc",
        <T extends ChartEventNames>(
            name: T,
            chartManager: ChartManager,
            value: ChartMsgValue<T>,
        ) => {
            if (name === "view_loaded") {
                //value type should be inferred as ViewId, but it's not?
                //(in contrast to below, where `if (datatype === 'popout')` works as expected)
                console.log(
                    "view_loaded",
                    value,
                    Object.keys(chartManager.charts),
                );
                if (!initialized) {
                    setupVuplex();
                    initialized = true;
                }
            }
        },
    );
    //temporarily disable websocket connection
    const socket = { on: (s: any, f: any) => {} }; //io(url);

    function sendMessage(msg: MDVMessage) {
        // socket.emit("message", msg);
        if (window.vuplex) window.vuplex.postMessage(msg);
    }

    const originalPopOutChart = cm._popOutChart;
    cm._popOutChart = (chart) => {
        originalPopOutChart.apply(cm, [chart]);
        sendMessage({ type: "popout", chartID: chart.config.id });
    };

    // socket.connect();
    function popout(chartID: string) {
        const chart = cm.charts[chartID];
        if (!chart) {
            console.error(`Chart ${chartID} not found`);
            console.log(
                "current chart keys: ",
                Object.keys(cm.charts).join(", "),
            );
            sendMessage({
                type: "error",
                message: `Popout failed: chart ${chartID} not found... trying again...`,
            });
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
        dataModel.addListener("IPC", () => {
            sendMessage({
                type: "filter",
                dataSource: ds.name,
                indices: dataModel.data,
            });
        });
    }

    function setupVuplex() {
        console.log("setupVuplex...");
        function addMessageListener() {
            if (!window.vuplex) {
                console.log("vuplex not found, cannot add listener (this was only relevant to unity embedded prototype)");
                return;
            }
            window.vuplex.postMessage({ type: "vuplex_ready" });
            window.vuplex.addEventListener("message", (e) => {
                let msg = e.data;
                if (typeof msg === "string")
                    msg = JSON.parse(msg) as MDVMessage;
                console.log("vuplex message", JSON.stringify(msg, null, 2));
                // TOOD: Types, Zod?
                try {
                    const data = msg; //JSON.parse(msg) as MDVMessage;
                    if (data.type === "popout") {
                        //TypeScript can infer from data.type === "popout" that data.chartID is a string based on what we've told it.
                        //(we don't *really* know that because we used 'as MDVMessage' without validating the data...
                        //so we could get a runtime error if the data is malformed, hence the comment above about Zod,
                        //but at least we get nice editor behavior and are explicit about the type we expect)
                        console.log("chartID: ", data.chartID);
                        popout(data.chartID);
                    }
                } catch (e) {
                    console.error("Error parsing vuplex message", e);
                    sendMessage({
                        type: "error",
                        message: `Error parsing vuplex message: ${e}`,
                    });
                }
            });
        }
        if (window.vuplex) {
            console.log("vuplex already exists, adding listener");
            addMessageListener();
        } else {
            console.log("vuplex doesn't exist (yet?), waiting for vuplexready");
            window.addEventListener("vuplexready", () => {
                console.log("vuplexready event received");
                addMessageListener();
            });
        }
    }
}
