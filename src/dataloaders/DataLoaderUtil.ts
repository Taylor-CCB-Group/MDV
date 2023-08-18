import type { Datasource } from "../modules/static_index";
import { getArrayBufferDataLoader, getLocalCompressedBinaryDataLoader } from "./DataLoaders";

export function getDataLoader(isStaticFolder: boolean, datasources: Datasource[], views: any, url: string) {
    const root = url.endsWith("/") ? url.substring(0, url.length-1) : url;

    return isStaticFolder ? {
        function: getLocalCompressedBinaryDataLoader(datasources, root),
        viewLoader: async (view: string) => views[view],
        rowDataLoader: loadRowDataStatic,
        binaryDataLoader: loadBinaryDataStatic
    } : {
        function: getArrayBufferDataLoader("/get_data"),
        viewLoader: getView,
        rowDataLoader: loadRowData,
        binaryDataLoader: loadBinaryData
    }
    //PJT - want to clarify /binarydata/ & /rowdata/ folders, and use of .b vs .gz

    async function loadRowDataStatic(datasource: string, index: string) {
        const resp = await fetch(`${root}/rowdata/${datasource}/${index}.json`);
        if (resp.status != 200) {
            return null
        }
        return await resp.json()
    }
    //load arbritray data
    async function loadBinaryDataStatic(datasource: string, name: string) {
        const resp = await fetch(`${root}/binarydata/${datasource}/${name}.b`);
        return await resp.arrayBuffer();
    }
}


//loads unstructured data for each row
async function loadRowData(datasource: string, index: string) {
    return await getPostData("/get_row_data", { datasource, index })
}
//load view from API
async function getView(view: string) {
    return await getPostData("/get_view", { view })
}

async function loadBinaryData(datasource: string, name: string) {
    return await getPostData("/get_binary_data", { datasource, name }, "arraybuffer");
}
//send json args and return json/array buffer response
export async function getPostData(url: string, args, return_type = "json") {
    const resp = await fetch(url, {
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
