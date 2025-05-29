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
type PingMessage = {
    type: "ping";
    message?: string;
}
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

type MDVMessage = PopoutMessage | FilterMessage | ErrorMessage | PingMessage;

/** initial experimental IPC support. As of this writing, will attempt to
 * - connect a 'vuplex' PostMessage interface for embedding in a VR environment in Unity (deprecated/no plans to use)
 * - connect a websocket to a socket.io server that will forward messages to other clients
 *    - started to use this for ai chatbot, started using REST but may move back
 *    - could be useful for other things, like syncing views, or receiving updates from a server
 *      e.g. about new data (e.g. from some kind of data pipeline), or about new views to load.
 */
export default async function connectIPC(cm: ChartManager) {
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
                    // setupVuplex();
                    initialized = true;
                }
            }
        },
    );
    // const socket = { on: (s: any, f: any) => {}, emit: (...args: any[])=>{} }; //io(url);
    // we should configure this with appropriate auth if needed, a namespace, etc.
    // const url = `${location.origin}/${location.pathname}`;
    // should use a better way of joining paths
    const socket = io(undefined, {path: `${cm.config.mdv_api_root}socket.io`});

    function sendMessage(msg: MDVMessage) {
        socket.emit("message", msg);
        if (window.vuplex) window.vuplex.postMessage(msg);
    }

    const originalPopOutChart = cm._popOutChart;
    cm._popOutChart = (chart) => {
        originalPopOutChart.apply(cm, [chart]);
        sendMessage({ type: "popout", chartID: chart.config.id });
    };

    await new Promise<void>((resolve, reject) => {
        socket.connect();
        socket.on("connect", () => {
            console.log("socket connected");
            sendMessage({type: 'ping'});
            // this is probably not the right design - what about disconnects?
            resolve();
        });
    });
    /** the first prototype feature built with this was making a window popout in response
     * to a message, so that in a VR environment in Unity we could see the chart in a separate window.
     */
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
    socket.on("message", async (data: MDVMessage) => {
        if (data.type === "popout") {
            const { chartID } = data;
            popout(chartID);
        }
        if (data.type === "ping") {
            console.log("ping received");
            console.log(data);
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
    return { socket, sendMessage };
}
