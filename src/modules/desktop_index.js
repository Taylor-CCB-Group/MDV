import "./all_css";
import ChartManager from "../charts/ChartManager.js";
import {
    getArrayBufferDataLoader,
    getLocalCompressedBinaryDataLoader,
    decompressData,
} from "../dataloaders/DataLoaders";
import { setProjectRoot } from "../dataloaders/DataLoaderUtil";
import { fetchJsonConfig } from "../dataloaders/DataLoaderUtil";

let route = "";

/**
 * @param {string=} routeFromTemplate - if undefined, it is assumed that the data is static
 * (i.e. the result of `convert_to_static_page()` in python), and located in the root of the current server.
 * - if it is the empty string, it will use the single-project API (`/get_data` etc).
 * - otherwise, it should be in the form of `/project/<project_id>`, which will be inserted by the Flask template
 * and use multi-project API (`/project/<project_id>/get_data` etc).
 */
function _mdvInit(routeFromTemplate) {
    const staticFolder = routeFromTemplate === undefined;
    if (!staticFolder) {
        route = routeFromTemplate;
        if (route && !route.startsWith("/project/")) {
            throw new Error(
                "routeFromTemplate must be undefined, the empty string, or in the form of '/project/<project_id>'",
            );
        }
    }
    setProjectRoot(route);
    //get the configs for MDV
    getConfigs(staticFolder).then((resp) => {
        const config = resp.state;
        
        //is view in the URL
        const urlParams = new URLSearchParams(window.location.search);
        const view = urlParams.get("view");
        if (config.all_views && view && config.all_views.indexOf(view) !== -1) {
            config["initial_view"] = view;
        }
        //data loaders depend on whether data is static or retrieved via API
        const dataLoader = {
            function: staticFolder
                ? getLocalCompressedBinaryDataLoader(resp.datasources, ".")
                : getArrayBufferDataLoader(`${route}/get_data`),
            viewLoader: staticFolder
                ? async (view) => resp.views[view]
                : getView,
            rowDataLoader: staticFolder ? loadRowDataStatic : loadRowData,
            binaryDataLoader: staticFolder
                ? loadBinaryDataStatic
                : loadBinaryData,
        };
        //listen to events
        const listener = (type, cm, data) => {
            switch (type) {
                case "state_saved":
                    //maybe consider rewriting the base URL here...
                    getData(`${route}/save_state`, data).then((resp) => {
                        if (resp.success) {
                            cm.createInfoAlert("Data Saved", {
                                duration: 2000,
                            });
                            cm.setAllColumnsClean();
                        } else {
                            cm.createInfoAlert("UnableToSaveData", {
                                duration: 3000,
                                type: "danger",
                            });
                        }
                    });
                    break;
                case "view_loaded":
                    changeURLParam("view", cm.viewManager.current_view);
                    break;
            }
        };
        //create the app
        console.log(config);
        new ChartManager(
            "holder",
            resp.datasources,
            dataLoader,
            config,
            listener,
        );
    });
}
//only method required
window._mdvInit = _mdvInit;

//loads unstructured data for each row
async function loadRowData(datasource, index) {
    return await getData(`${route}/get_row_data`, { datasource, index });
}
async function loadRowDataStatic(datasource, index) {
    const resp = await fetch(`${route}/rowdata/${datasource}/${index}.json`);
    if (resp.status !== 200) {
        return null;
    }
    return await resp.json();
}
//load view from API
async function getView(view) {
    return await getData(`${route}/get_view`, { view: view });
}

//load arbitrary data
async function loadBinaryDataStatic(datasource, name) {
    const resp = await fetch(
        `${route}/binarydata/${datasource}/${name}.gz, {responseType: "arraybuffer"}`,
    );
    const b = await resp.arrayBuffer();
    return await decompressData(b);
}
async function loadBinaryData(datasource, name) {
    const b = await getData(
        `${route}/get_binary_data`,
        { datasource, name },
        "arraybuffer",
    );
    return await decompressData(b);
}

/**
 * get the configs whether from a folder or via a remote API
 * @param {boolean} isStaticFolder - true if configs are in a folder, false for remote API
 * @returns {object} - configs (state, datasources, views)
 */
async function getConfigs(isStaticFolder) {
    let configs = {};
    const fetch = async (url) => fetchJsonConfig(url, route);
    if (isStaticFolder) {
        let resp = await fetch(`${route}/datasources.json`);
        configs.datasources = await resp.json();
        resp = await fetch(`${route}/state.json`);
        configs.state = await resp.json();
        resp = await fetch(`${route}/views.json`);
        configs.views = await resp.json();
    } else {
        configs = await getData(`${route}/get_configs`);
    }
    return configs;
}
/** send json args and return json/array buffer response */
async function getData(url, args, return_type = "json") {
    const resp = await fetch(url, {
        method: "POST",
        body: JSON.stringify(args),
        headers: {
            Accept: "application/json,text/plain,*/*",
            "Content-Type": "application/json",
        },
    });
    if (return_type === "json") {
        console.log("returning json, rewriting base_url first...");
        const original = await resp.json();
        console.log("original", original);
        return original; //rewriteBaseUrlRecursive(original, route);
    }

    return await resp.arrayBuffer();
}

//changes or adds a param to the browser address bar
export function changeURLParam(param, value) {
    const url = new URL(window.location);
    url.searchParams.has(param)
        ? url.searchParams.set(param, value)
        : url.searchParams.append(param, value);
    url.search = url.searchParams;
    history.pushState({}, null, url);
}
