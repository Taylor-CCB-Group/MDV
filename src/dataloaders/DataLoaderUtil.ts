import type { Datasource } from "../modules/static_index";
import { getArrayBufferDataLoader, getLocalCompressedBinaryDataLoader } from "./DataLoaders";

export function getDataLoader(isStaticFolder: boolean, datasources: Datasource[], views: any, url: string) {
    const root = url.endsWith("/") ? url.substring(0, url.length-1) : url;

    function rewriteBaseUrlRecursive(config) {
        if (root === "./") return;
        if (Array.isArray(config)) {
            for (const item of config) {
                rewriteBaseUrlRecursive(item);
            }
            return;
        }
        for (const key in config) {
            if (key === "base_url") {
                config[key] = config[key].replace("./", `${root}/`);
            } else if (typeof config[key] === "object") {
                rewriteBaseUrlRecursive(config[key]);
            }
        }
    }
    
    async function fetchAndPatchJSON(url: string) {
        let resp = await fetch(url);
        const config = await resp.json();
        rewriteBaseUrlRecursive(config);
        return config;
    }


    return isStaticFolder ? {
        function: getLocalCompressedBinaryDataLoader(datasources, root),
        viewLoader: async (view: string) => views[view],
        rowDataLoader: async (dsName: string) => await fetchAndPatchJSON(`${root}/${dsName}.json`),
        binaryDataLoader: loadBinaryDataStatic
    } : {
        function: getArrayBufferDataLoader("/get_data"),
        viewLoader: getView,
        rowDataLoader: loadRowData,
        binaryDataLoader: loadBinaryData
    }
}


//loads unstructured data for each row
async function loadRowData(datasource, index) {
    return await getPostData("/get_row_data", { datasource, index })
}
async function loadRowDataStatic(datasource, index) {
    const resp = await fetch(`./rowdata/${datasource}/${index}.json`);
    if (resp.status != 200) {
        return null
    }
    return await resp.json()
}
//load view from API
async function getView(view) {
    return await getPostData("/get_view", { view: view })
}

//load arbritray data
async function loadBinaryDataStatic(datasource, name) {
    const resp = await fetch(`./binarydata/${datasource}/${name}.b`);
    return await resp.arrayBuffer();
}
async function loadBinaryData(datasource, name) {
    return await getPostData("/get_binary_data", { datasource, name }, "arraybuffer");
}
//send json args and return json/array buffer response
async function getPostData(url: string, args, return_type = "json") {
    const resp = await fetch(url,
        {
            method: "POST",
            body: JSON.stringify(args),
            headers: {
                "Accept": "application/json,text/plain,*/*",
                "Content-Type": "application/json"
            }
        });
    if (return_type === "json") {
        return await resp.json();
    }
    else {
        return await resp.arrayBuffer();
    }
}
