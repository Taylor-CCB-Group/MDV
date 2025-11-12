/// <reference types="vite/client" />
import "./wdyr";
import "@slickgrid-universal/common/dist/styles/css/slickgrid-theme-material.css";
import "./all_css";
import HmrHack from "../react/HmrHack";
HmrHack();
import ChartManager from "../charts/ChartManager.js";

import {
    fetchJsonConfig,
    getDataLoader,
    setProjectRoot,
} from "../dataloaders/DataLoaderUtil";
import { changeURLParam } from "./desktop_index";
import BaseChart from "../charts/BaseChart";
import type { DataSource } from "@/charts/charts";
import { getProjectName } from "./ProjectContext";
import { createMdvPortal } from "@/react/react_utils";
import ProjectStateHandlerWrapper from "@/react/ProjectStateHandler";
import type { Root } from "react-dom/client";

// see also basic_index.js for some global mdv stuff... only interested in chartManager for now.
declare global {
    interface Window {
        mdv: {
            ChartManager: typeof ChartManager;
            chartManager: ChartManager;
            chartTypes?: any;
            debugChart?: any;
        };
    }
}
//@ts-expect-error maybe we'll be less hacky about this later
window.mdv = {
    // does this not want to be a singleton - rather than referring to the class?
    // not sure what I'd use the class for - would generally import in a module
    ChartManager,
    chartTypes: BaseChart.types,
};

document.addEventListener("DOMContentLoaded", async () => {
    try {
        await loadData();
    } catch (e) {
        // alert(`Error loading data: ${e}`)
        document.head.title = "MDV - Error loading data"; // not actually appearing, hey-ho.
        document.body.style.textAlign = "center";
        document.body.innerHTML = `<h1>Error loading data</h1><p>${e}</p>`;
    }
});

/// --- this section is a bit of a mess, but it works for now ---
//(copilot suggesting "this is a bit of a mess" is a bit rude, but in this case it's an understatement)
// todo: change so that *only* ?dir=... is used to determine the root?
// const flaskURL = window.location.pathname;
const { origin, pathname } = window.location;
const flaskURL = origin + pathname; //new URL(href).origin;//href.substring(href.indexOf("/project"));
// if URLSearchParams has a 'dir' parameter, use that as the data directory.
const urlParams = new URLSearchParams(window.location.search);
// if we're in a popout window, ignore the dir parameter and don't load data
const isPopout = urlParams.get("popout") === "true";
// if there is no dir parameter, use the flaskURL to proxy requests to the python server
const dir = urlParams.get("dir") || (isPopout ? "" : flaskURL);
function getRoot(dir: string) {
    // const url = new URL(dir);
    // return url.origin + url.pathname;
    return dir.endsWith("/") ? dir.substring(0, dir.length - 1) : dir;
}
const root = getRoot(dir);
//hack to load data from local API... TODO XXX make less hacky, rather than more...
//this was sort-of-working as of this writing for `MDVProject.serve()`, as long as the default port 5050 is used...
//need to revisit and actually come up with a proper design.
// - this means some online projects are currently broken, I should fix that. <--
// as of this writing, they work with a `?static` parameter, but that's not a good solution.
const staticFolder = urlParams.get("static") !== null; //!dir.startsWith("/project") && !(window.location.port === "5050") && !dir.endsWith("5050");
const project_id = dir.split("/").pop();
// State handler container and root
let stateHandlerContainer: HTMLElement | null = null;
let stateHandlerRoot: Root | null = null;

// getting the project name by passing project id
getProjectName(Number(project_id)).then((project_name) => {
    document.title = `MDV - ${project_name}`
});
/// --- end of messy section ---

// set title of page to the data directory
if (isPopout) document.title = "MDV popout";

async function loadData() {
    // setupDebug();
    if (isPopout) return;
    setProjectRoot(root);
    // move more of this into DataLoaderUtil (which might get renamed)?
    const datasources = (await fetchJsonConfig(
        `${root}/datasources.json`,
        root,
        true,
    )) as DataSource[];
    const config = await fetchJsonConfig(`${root}/state.json`, root, true);
    config.popouturl = undefined;
    // todo: check if this is correct
    const permission = config?.permission === "edit" || config?.permission === "owner";
    const views = await fetchJsonConfig(`${root}/views.json`, root, true);
    //is view in the URL
    const view = urlParams.get("view");
    if (config.all_views && view && config.all_views.indexOf(view) !== -1) {
        config["initial_view"] = view;
    }
    const isStatic = staticFolder || config.static || false;

    const dataLoader = getDataLoader(isStatic, datasources, views, dir);

    const listener = async (type: string, cm: ChartManager, data: any) => {
        if (type === "state_saved") {
            // Unmount the existing root
            if (stateHandlerRoot) {
                stateHandlerRoot.unmount();
            }

            // Create a new container if it's already not created
            if (!stateHandlerContainer) {
                stateHandlerContainer = document.createElement('div');
                stateHandlerContainer.id = "mdv_state_handler";
                document.body.appendChild(stateHandlerContainer);
            }
            stateHandlerRoot = createMdvPortal(ProjectStateHandlerWrapper({root, data, staticFolder, permission}), stateHandlerContainer);
        }
        if (type === "view_loaded") {
            changeURLParam("view", cm.viewManager.current_view);
        }
    };
    // this will assign itself to `window.mdv.chartManager` in constructor
    const cm = new ChartManager(
        "holder",
        datasources,
        dataLoader,
        config,
        listener as any, //jsdoc 🙄
    );
}
