import type { DataSource } from "@/charts/charts";
import { fetchJsonConfig, getDataLoader, setProjectRoot } from "@/dataloaders/DataLoaderUtil";
import { getApiRootFromDir, getProjectDirFromLocation } from "@/utils/mdvRouting";

export type ProjectBootstrapContext = {
    urlParams: URLSearchParams;
    isPopout: boolean;
    dir: string;
    root: string;
    staticFolder: boolean;
    projectId: string;
    apiRoot: string;
};

export type LoadedProjectRuntime = {
    datasources: DataSource[];
    config: any;
    views: any;
    dataLoader: ReturnType<typeof getDataLoader>;
    permission: boolean;
    staticFolder: boolean;
    root: string;
    dir: string;
};

function normalizeRoot(dir: string) {
    return dir.endsWith("/") ? dir.slice(0, -1) : dir;
}

export function getProjectBootstrapContext(): ProjectBootstrapContext {
    const urlParams = new URLSearchParams(window.location.search);
    const isPopout = urlParams.get("popout") === "true";
    const dir = isPopout ? "" : getProjectDirFromLocation();
    const root = normalizeRoot(dir);
    const staticFolder = urlParams.get("static") !== null;
    const projectId = dir.split("/").pop() ?? "";
    const apiRoot = getApiRootFromDir(root);

    return { urlParams, isPopout, dir, root, staticFolder, projectId, apiRoot };
}

export async function loadProjectRuntime(
    context: ProjectBootstrapContext,
): Promise<LoadedProjectRuntime> {
    const { root, dir, staticFolder, urlParams } = context;
    setProjectRoot(root);

    const [datasources, config, views] = await Promise.all([
        fetchJsonConfig<DataSource[]>(`${root}/datasources.json`, root, false),
        fetchJsonConfig(`${root}/state.json`, root, false),
        fetchJsonConfig(`${root}/views.json`, root, false),
    ]);

    config.popouturl = undefined;
    const permission = config?.permission === "edit";
    const view = urlParams.get("view");
    if (config.all_views && view && config.all_views.indexOf(view) !== -1) {
        config.initial_view = view;
    }

    const isStatic = staticFolder || config.static || false;
    const dataLoader = getDataLoader(isStatic, datasources, views, dir);

    return {
        datasources,
        config,
        views,
        dataLoader,
        permission,
        staticFolder: isStatic,
        root,
        dir,
    };
}
