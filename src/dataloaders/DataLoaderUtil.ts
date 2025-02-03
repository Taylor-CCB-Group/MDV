import { createEl } from "@/utilities/ElementsTyped";
import {
    getArrayBufferDataLoader,
    getLocalCompressedBinaryDataLoader,
} from "./DataLoaders";
import type { DataSource, DataColumn, DataType, LoadedDataColumn } from "@/charts/charts";
import { isColumnLoaded } from "@/lib/columnTypeHelpers";

let projectRoot = "";
export async function fetchJsonConfig(url: string, root: string) {
    try {
        const resp = await fetch(url);
        const config = await resp.json();
        //rewriteBaseUrlRecursive(config, root); //removed.
        return config;
    } catch (e) {
        console.error(`Error fetching ${url}: ${e}`);
        return { error: e };
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

export function getDataLoader(
    isStaticFolder: boolean,
    datasources: DataSource[],
    views: any,
    url: string,
) {
    const root = url.endsWith("/") ? url.substring(0, url.length - 1) : url;

    return isStaticFolder
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
    //PJT - want to clarify /binarydata/ & /rowdata/ folders, and use of .b vs .gz

    async function loadRowDataStatic(datasource: string, index: string) {
        const resp = await fetch(`${root}/rowdata/${datasource}/${index}.json`);
        if (resp.status !== 200) {
            return null;
        }
        return await resp.json();
    }
    //load arbritray data
    async function loadBinaryDataStatic(datasource: string, name: string) {
        const resp = await fetch(`${root}/binarydata/${datasource}/${name}.b`);
        return await resp.arrayBuffer();
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
        //todo nicer dialog
        const dialog = createEl("dialog", { title: "Error loading view" });
        dialog.innerText = e.message;
        dialog.addEventListener("click", () => dialog.close());
        document.body.appendChild(dialog);
        dialog.showModal();
        console.error(e);
        // return {initialCharts: {}};
        return undefined;
    }
}

async function loadBinaryData(datasource: string, name: string) {
    return await getPostData(
        `${projectRoot}/get_binary_data`,
        { datasource, name },
        "arraybuffer",
    );
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
