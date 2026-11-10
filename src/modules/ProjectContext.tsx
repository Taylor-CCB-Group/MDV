import axios from "axios";
import { createContext, type PropsWithChildren, useContext, useEffect, useState } from "react";
import useBuildInfo, { type BuildInfo } from "@/catalog/hooks/useBuildInfo";
import type ChartManager from "@/charts/ChartManager";
import { buildApiUrl, getApiRootFromDir, getProjectDirFromLocation, normalizeApiRoot } from "@/utils/mdvRouting";

export type ProjectInfoBase = {
    root: string;
    staticFolder: boolean;
    chartManager: ChartManager;
    projectName: string;
    projectApiRoute: string;
    mainApiRoute: string;
};

export type ProjectInfo = ProjectInfoBase & {
    buildInfo: BuildInfo;
};

/** One project, as `GET /projects/<id>` returns it. */
export type Project = {
    id: number;
    name: string;
    lastModified: string;
    thumbnail?: string | null;
    readme?: string | null;
    owner: string[];
    writable: boolean;
    /** Null when auth is off. */
    permissions: { can_read: boolean; can_write: boolean; is_owner: boolean } | null;
};

/**
 * Get project information from the URL and URL parameters.
 *
 * This can be used to derive the initial state of the project, or called by non-React code to get similar information.
 */
function getProjectInfoBaseFromChartManager(chartManagerOverride?: ChartManager): ProjectInfoBase {
    // Derive initial state from window.location and URL parameters
    const urlParams = new URLSearchParams(window.location.search);
    const dir = getProjectDirFromLocation();

    // let's adjust how we reason about root for main API and project API
    // also consider 'static' - when we output state.json for that
    const getRoot = (dir: string) => {
        return dir.endsWith("/") ? dir.substring(0, dir.length - 1) : dir;
    };
    const root = getRoot(dir);
    // todo - get these from e.g. state.json instead?
    const staticFolder = urlParams.get("static") !== null;
    // this is really projectID in the current implementation - would be good to have name as well (and be able to change it)
    const projectName = dir.split("/").pop() || ""; //todo - check logic for default project name
    const chartManager = chartManagerOverride ?? window.mdv?.chartManager;
    if (!chartManager) {
        throw new Error("ChartManager is not initialised.");
    }
    //! might need to be a bit careful if we end up using this
    const mainApiRoute = normalizeApiRoot(chartManager.config.mdv_api_root || getApiRootFromDir(root));
    // const root = mainApiRoute; // todo establish that we have an actual consistent logic for this
    //! only applies when running with the associated project API routing...
    const projectApiRoute = buildApiUrl(`project/${projectName}`, mainApiRoute).replace(/\/$/, "") + "/";
    return {
        root,
        staticFolder,
        projectName,
        chartManager,
        mainApiRoute,
        projectApiRoute,
    };
}

export function getProjectInfoBase(): ProjectInfoBase {
    return getProjectInfoBaseFromChartManager();
}

export function useProjectInfo(): ProjectInfo {
    const { buildInfo } = useBuildInfo();
    return {
        ...getProjectInfoBase(),
        buildInfo,
    };
}

/** Fetch one project from `GET /projects/<id>`, or null if the server doesn't have it. */
export async function getProject(projectId: number, apiRoot?: string): Promise<Project | null> {
    try {
        const res = await axios.get<Project>(buildApiUrl(`projects/${projectId}`, apiRoot));
        return res.data;
    } catch (error) {
        console.error("Failed to fetch project:", error);
        return null;
    }
}

/** The project this page is showing, or null until it loads or if it can't be found. */
export function useCurrentProject(): Project | null {
    const { projectName: projectId, root, mainApiRoute } = useProject();
    const apiRoot = mainApiRoute === "/" ? getApiRootFromDir(root) : mainApiRoute;
    const [project, setProject] = useState<Project | null>(null);
    useEffect(() => {
        let cancelled = false;
        void getProject(Number(projectId), apiRoot).then((result) => {
            if (!cancelled) setProject(result);
        });
        return () => {
            cancelled = true;
        };
    }, [projectId, apiRoot]);
    return project;
}

const ProjectContext = createContext<ProjectInfo>(null as any);

export const ProjectProvider = ({ children, chartManager }: PropsWithChildren<{ chartManager?: ChartManager }>) => {
    const { buildInfo } = useBuildInfo();
    const base = getProjectInfoBaseFromChartManager(chartManager);

    return <ProjectContext.Provider value={{ ...base, buildInfo }}>{children}</ProjectContext.Provider>;
};

export const useProject = () => useContext(ProjectContext);
