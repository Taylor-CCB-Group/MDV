import type ChartManager from "@/charts/ChartManager";
import React, {
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

    const getRoot = (dir: string) => {
        return dir.endsWith("/") ? dir.substring(0, dir.length - 1) : dir;
    };
    const root = getRoot(dir);
    // todo - get these from e.g. state.json instead?
    const staticFolder = urlParams.get("static") !== null;
    // this is really projectID in the current implementation - would be good to have name as well (and be able to change it)
    const projectName = dir.split("/").pop() || ""; //todo - check logic for default project name
    const { chartManager } = window.mdv;
    return {
        root,
        staticFolder,
        projectName,
        chartManager,
    };   
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
