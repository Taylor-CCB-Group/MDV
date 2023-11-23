import type { Datasource } from "../modules/static_index";
import { getArrayBufferDataLoader, getLocalCompressedBinaryDataLoader } from "./DataLoaders";

/**
 * mutate the config object to substitute the project root into all base_url entries where relevant.
 * returns the config object.
 */
export function rewriteBaseUrlRecursive(config, root: string) {
    console.warn("rewriteBaseUrlRecursive - no-op: use getProjectURL() instead / subject to review.", config, root)
    return config;
    if (Array.isArray(config)) {
        for (const item of config) {
            rewriteBaseUrlRecursive(item, root);
        }
        return config;
    }
    for (const key in config) {
        if (key === "base_url") {
            // when a project has been saved, it can have a base_url that doesn't start with "./", but also doesn't match our root...
            // e.g. 'http://localhost:5050/images/...'.
            // we may need a better spec for thinking about this.
            // one hacky solution for now *may* be to make sure we only use these modified base_urls internally, and never save them.
            // which we could possibly do by doing the oposite of this function when saving...
            // that seems potentially fragile though.
            // there aren't many places that use base_url, so we could just change them all to use root instead.
            // want to be a bit careful in case there could at some point be multiple projects in the same global js context.
            const newUrl = config[key].replace("./", `${root}/`);
            // console.log("rewriting base url", config[key], newUrl);
            // config[key] = config[key].replace("./", newUrl);
            // rather than replacing the base_url, we could add a new url_in_project key...
            // which will get saved along with the base_url, but will be used internally instead.
            config['url_in_project'] = newUrl;
        } else if (typeof config[key] === "object") {
            rewriteBaseUrlRecursive(config[key], root);
        }
    }
    return config;
}
let projectRoot = "";
export async function fetchAndPatchJSON(url: string, root: string) {
    let resp = await fetch(url);
    const config = await resp.json();
    rewriteBaseUrlRecursive(config, root);
    return config;
}

export function setProjectRoot(root: string) {
    projectRoot = root;
}

//https://stackoverflow.com/questions/29855098/is-there-a-built-in-javascript-function-similar-to-os-path-join
export function buildPath(...args: string[]) {
    return args.map((part, i) => {
        if (i === 0) {
            return part.trim().replace(/[\/]*$/g, '')
        } else {
            return part.trim().replace(/(^[\/]*|[\/]*$)/g, '')
        }
    }).filter(x => x.length).join('/')
}

/**
 * Given a url that is relative to the project root, return a more API-appropriate URL.
 */
export function getProjectURL(url: string) {
    if (url.startsWith(projectRoot)) return url;
    return buildPath(projectRoot, url).replace("./", "/") + '/';
}


export function getDataLoader(isStaticFolder: boolean, datasources: Datasource[], views: any, url: string) {
    const root = url.endsWith("/") ? url.substring(0, url.length-1) : url;

    return isStaticFolder ? {
        function: getLocalCompressedBinaryDataLoader(datasources, root),
        viewLoader: async (view: string) => views[view],
        rowDataLoader: loadRowDataStatic,
        binaryDataLoader: loadBinaryDataStatic
    } : {
        function: getArrayBufferDataLoader(root+"/get_data"),
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
    return await getPostData(projectRoot+"/get_row_data", { datasource, index })
}
//load view from API
async function getView(view: string) {
    return await getPostData(projectRoot+"/get_view", { view })
}

async function loadBinaryData(datasource: string, name: string) {
    return await getPostData(projectRoot+"/get_binary_data", { datasource, name }, "arraybuffer");
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
