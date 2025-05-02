import type ChartManager from "@/charts/ChartManager";
import axios from "axios";
import {
    type PropsWithChildren,
    createContext,
    useContext,
    useState,
} from "react";

export type BuildInfo = {
    commitDate: string;
    branchName: string;
    commitHash: string;
    lastCommitMessage: string;
    buildDate: string;
    dirty: boolean;
}

export type ProjectInfo = {
    root: string;
    staticFolder: boolean;
    chartManager: ChartManager;
    projectName: string;
    // buildInfo: BuildInfo;
};

export type Project = {
    id: number;
    lastModified: string;
    name: string;
};

/**
 * Returns information about the current build, particularly git revision.
 * 
 * Note, hardcoded for vite, this will need to be updated if we switch to a different bundler.
 */
export function getBuildInfo(): BuildInfo {
    const commitDate = import.meta.env.VITE_GIT_COMMIT_DATE || "";
    const branchName = import.meta.env.VITE_GIT_BRANCH_NAME || "";
    const commitHash = import.meta.env.VITE_GIT_COMMIT_HASH || "";
    const lastCommitMessage = import.meta.env.VITE_GIT_LAST_COMMIT_MESSAGE || "";
    const buildDate = import.meta.env.VITE_BUILD_DATE || "";
    const dirty = import.meta.env.VITE_GIT_DIRTY === "dirty";
    return {
        commitDate,
        branchName,
        commitHash,
        lastCommitMessage,
        buildDate,
        dirty,
    };
}

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
        // buildInfo: getBuildInfo(),
    };   
}

export async function getProjectName(projectId: number) {
    try {
        const res = await axios.get("/projects");
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
