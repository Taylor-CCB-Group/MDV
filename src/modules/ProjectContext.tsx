import type ChartManager from "@/charts/ChartManager";
import React, {
    type PropsWithChildren,
    createContext,
    useContext,
    useState,
} from "react";

type ProjectInfo = {
    root: string;
    staticFolder: boolean;
    chartManager: ChartManager;
    projectName: string;
};

const ProjectContext = createContext<ProjectInfo>(null as any);

export const ProjectProvider = ({ children }: PropsWithChildren) => {
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
    const projectName = dir.split("/").pop() || ""; //todo - check logic for default project name
    const { chartManager } = window.mdv;

    const [projectConfig, setProjectConfig] = useState<ProjectInfo>({
        // not used? - test / check <<<---
        // I think we determined that `root` was the useful thing to expose
        // in a way that hides most of the sketchiness of the current setup.
        // flaskURL,
        // dir,
        root,
        staticFolder,
        projectName,
        chartManager,
    });

    return (
        <ProjectContext.Provider value={{ ...projectConfig }}>
            {children}
        </ProjectContext.Provider>
    );
};

export const useProject = () => useContext(ProjectContext);
