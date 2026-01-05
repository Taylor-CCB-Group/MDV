import {
    getArrayBufferDataLoader,
    getLocalCompressedBinaryDataLoader,
} from "./DataLoaders";
import type { DataSource, DataColumn, DataType, LoadedDataColumn } from "@/charts/charts";
import { isColumnLoaded } from "@/lib/columnTypeHelpers";
import { createMdvPortal } from "@/react/react_utils";
import ErrorComponentReactWrapper from "@/react/components/ErrorComponentReactWrapper";
import { decompressData } from "./DataLoaders";

let projectRoot = "";
/**
 * This could potentially also have some more awareness of what type of json it is fetching, and e.g. do some zod validation
 */
export async function fetchJsonConfig(url: string, root: string, createErrorComponent = false)  {
    try {
        const resp = await fetch(url);
        const config = await resp.json();
        if (!resp.ok) {
            throw new Error(JSON.stringify(config, null, 2));
        }
        //rewriteBaseUrlRecursive(config, root); //removed.
        return config;
    } catch (error: any) {
        if (createErrorComponent) {
            // Use a dedicated container instead of document.body to avoid React warnings
            let errorContainer = document.getElementById("mdv-error-container");
            if (!errorContainer) {
                errorContainer = document.createElement("div");
                errorContainer.id = "mdv-error-container";
                document.body.appendChild(errorContainer);
            }
            // Basic styling for the error container
            errorContainer.style.position = "fixed";
            errorContainer.style.top = "0";
            errorContainer.style.left = "0";
            errorContainer.style.width = "100%";
            errorContainer.style.height = "100%";
            errorContainer.style.display = "flex";
            errorContainer.style.justifyContent = "center";
            errorContainer.style.alignItems = "center";
            errorContainer.style.zIndex = "9999";
            errorContainer.style.backgroundColor = "rgba(0,0,0,0.1)";

            createMdvPortal(ErrorComponentReactWrapper({ error: {message: `Error fetching JSON '${url}'`}, extraMetaData: {message: `${error}`} }), errorContainer);
        }
        console.error(`Error fetching ${url}: ${error}`);
        throw error;
    }
}

export function setProjectRoot(root: string) {
    projectRoot = root;
}

//https://stackoverflow.com/questions/29855098/is-there-a-built-in-javascript-function-similar-to-os-path-join
export function buildPath(...args: string[]) {
    return args
        .map((part, i) => {
            if (i === 0) {
                return part.trim().replace(/[\/]*$/g, "");
            }
            return part.trim().replace(/(^[\/]*|[\/]*$)/g, "");
        })
        .filter((x) => x.length)
        .join("/");
}

/**
 * Given a url that is relative to the project root, return a more API-appropriate URL.
 *
 * ==This needs a rethink.==
 * Been meaning to change to `project=...` in the URL, especially now that `dir` gets rewritten
 * and ugli-fied.
 *
 * - where does this get applied? I missed `feature.js` - are there other places?
 * - there may be occasions where a project is setup with references to external resources
 */
export function getProjectURL(url: string, trailingSlash = true) {
    if (url.startsWith(projectRoot)) return url; //<<--- review this, test calling it multiple times
    const p =
        buildPath(projectRoot, url).replace("./", "/") +
        (trailingSlash ? "/" : "");
    // https://github.com/hms-dbmi/viv/issues/756
    if (p.startsWith("/")) return new URL(p, window.location.href).href;
    return p;
}
export type DataLoader = {
    function: (columns: DataColumn<any>[], dataSource: string, size: number) => Promise<ArrayBufferLike>,
    viewLoader: (view: string) => Promise<any>,
    rowDataLoader: (dataSource: string, index: string) => Promise<any>,
    binaryDataLoader: (dataSource: string, name: string) => Promise<ArrayBufferLike>,
}

export function getDataLoader(
    isStaticFolder: boolean,
    datasources: DataSource[],
    views: any,
    url: string,
) {
    const root = url.endsWith("/") ? url.substring(0, url.length - 1) : url;
    const loader = isStaticFolder
        ? {
              function: getLocalCompressedBinaryDataLoader(datasources, root),
              viewLoader: async (view: string) => views[view],
              rowDataLoader: loadRowDataStatic,
              binaryDataLoader: loadBinaryDataStatic,
          }
        : {
              function: getArrayBufferDataLoader(`${root}/get_data`),
              viewLoader: getView,
              rowDataLoader: loadRowData,
              binaryDataLoader: loadBinaryData,
          };
    return loader;
    //PJT - want to clarify /binarydata/ & /rowdata/ folders, and use of .b vs .gz

    async function loadRowDataStatic(datasource: string, index: string) {
        const resp = await fetch(`${root}/rowdata/${datasource}/${index}.json`);
        if (resp.status !== 200) {
            return null;
        }
        return await resp.json();
    }
    //load arbitrary data
    async function loadBinaryDataStatic(datasource: string, name: string) {
        const resp = await fetch(`${root}/binarydata/${datasource}/${name}.gz`);
        const buff =  await resp.arrayBuffer();
        return await decompressData(new Uint8Array(buff));
    }
}

//loads unstructured data for each row
async function loadRowData(datasource: string, index: string) {
    return await getPostData(`${projectRoot}/get_row_data`, {
        datasource,
        index,
    });
}
//load view from API
async function getView(view?: string) {
    if (!view) return undefined; //this seems to be a somewhat reasonable way to handle empty project with no views - but could do with more testing...
    try {
        return await getPostData(`${projectRoot}/get_view`, { view });
    } catch (e: any) {
        //todo add an error dialog here to display error
        console.error(e);
        // return {initialCharts: {}};
        return undefined;
    }
}

async function loadBinaryData(datasource: string, name: string) {
    const data =  await getPostData(
        `${projectRoot}/get_binary_data`,
        { datasource, name },
        "arraybuffer",
    );

    return await decompressData(data);
}
//send json args and return json/array buffer response
export async function getPostData(url: string, args: any, return_type = "json") {
    const resp = await fetch(url, {
        method: "POST",
        body: JSON.stringify(args),
        headers: {
            Accept: "application/json,text/plain,*/*",
            "Content-Type": "application/json",
        },
    });
    if (!resp.ok) {
        throw new Error(`Error fetching '${url}': ${resp.statusText}`);
    }
    if (return_type === "json") {
        return await resp.json();
    }

    return await resp.arrayBuffer();
}

/** Get a column with given name from the given dataSource, fetching from server if necessary.
 * Returns a Column object when the data is loaded.
 * 
 * todo - consider batching requests (add to queue, then load all at once in `setTimeout(...,0)`?)
 * better types - for `columnName` and return value generic(?)
 */
export async function loadColumn(datasourceName: string, columnName: string): Promise<LoadedDataColumn<DataType>> {
    return new Promise((resolve, reject) => {
        try {
            const ds = window.mdv.chartManager.getDataSource(datasourceName);
            const column = ds.columnIndex[columnName];
            if (!column) throw `no columnIndex['${columnName}'] record`;
            if (ds.columnsWithData.includes(columnName)) {
                if (!isColumnLoaded(column)) throw "assertion failed..."
                resolve(column); //hopefully this is trustworthy
            } else {
                window.mdv.chartManager.loadColumnSet(
                    [columnName],
                    datasourceName,
                    (failedColumns: DataColumn<any>[]) => {
                        if (!isColumnLoaded(column)) throw "broken promise for loading column"
                        if (failedColumns.length) {
                            // reject(new Error(`Failed to load column ${columnName}`));
                            console.error(`Failed to load column ${columnName}`);
                        }
                        return resolve(column);
                    }
                );
            }
        } catch (e) {
            reject(e);
        }
    });
}
