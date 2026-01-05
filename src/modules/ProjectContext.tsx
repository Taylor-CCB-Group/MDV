import useBuildInfo, {type BuildInfo} from "@/catalog/hooks/useBuildInfo";
import type ChartManager from "@/charts/ChartManager";
import axios from "axios";
import {
    type PropsWithChildren,
    createContext,
    useContext,
    useState,
} from "react";

export type ProjectInfo = {
    root: string;
    staticFolder: boolean;
    chartManager: ChartManager;
    projectName: string;
    buildInfo: BuildInfo;
    projectApiRoute: string;
    mainApiRoute: string;
};

export type Project = {
    id: number;
    lastModified: string;
    name: string;
};

/** 
 * Get project information from the URL and URL parameters.
 * 
 * This can be used to derive the initial state of the project, or called by non-React code to get similar information.
 */
export function getProjectInfo(): ProjectInfo {
    // Derive initial state from window.location and URL parameters
    const { origin, pathname } = window.location;
    const flaskURL = origin + pathname;
    const urlParams = new URLSearchParams(window.location.search);
    const dir = urlParams.get("dir") || flaskURL;

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
    const { chartManager } = window.mdv || {};
    //! might need to be a bit careful if we end up using this
    const mainApiRoute = chartManager?.config?.mdv_api_root || "/";
    // const root = mainApiRoute; // todo establish that we have an actual consistent logic for this
    //! only applies when running with the associated project API routing...
    const projectApiRoute = `${mainApiRoute}/project/${projectName}/`.replace('//', '/');
    const { buildInfo } = useBuildInfo();
    return {
        root,
        staticFolder,
        projectName,
        chartManager,
        buildInfo,
        mainApiRoute,
        projectApiRoute,
    };   
}

export async function getProjectName(projectId: number) {
    try {
        const res = await axios.get("../projects");
        const projectData: Project[] = res.data;
    
        const projectName = projectData?.find((project) => project.id === projectId)?.name || "unnamed_project";
        return projectName;
    } catch (error) {
        console.error("Failed to fetch project name:", error);
        return "unnamed_project";
    }
}

const ProjectContext = createContext<ProjectInfo>(null as any);

export const ProjectProvider = ({ children }: PropsWithChildren) => {
    // we might get away with just a hook for `useProject = () => useMemo(getProjectInfo, [])`... 
    // but maybe in future we'll make this more complex
    const [projectConfig, setProjectConfig] = useState<ProjectInfo>(getProjectInfo());

    return (
        <ProjectContext.Provider value={{ ...projectConfig }}>
            {children}
        </ProjectContext.Provider>
    );
};

export const useProject = () => useContext(ProjectContext);
