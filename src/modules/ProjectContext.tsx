import type ChartManager from '@/charts/ChartManager';
import React, { type PropsWithChildren, createContext, useContext, useState } from 'react';

type ProjectInfo = {
    root: string;
    staticFolder: boolean;
    chartManager: ChartManager;
}

const ProjectContext = createContext<ProjectInfo>(null);

export const ProjectProvider = ({ children }: PropsWithChildren) => {
    // Derive initial state from window.location and URL parameters
    const { origin, pathname } = window.location;
    const flaskURL = origin + pathname;
    const urlParams = new URLSearchParams(window.location.search);
    const dir = urlParams.get('dir') || flaskURL;

    const getRoot = (dir: string) => {
        return dir.endsWith("/") ? dir.substring(0, dir.length - 1) : dir;
    };
    const root = getRoot(dir);
    // todo - get these from e.g. state.json instead?
    const staticFolder = urlParams.get('static') !== null;
    const projectName = dir.split("/").pop();
    const { chartManager } = window.mdv;

    const [projectConfig, setProjectConfig] = useState({
        flaskURL,
        dir,
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
