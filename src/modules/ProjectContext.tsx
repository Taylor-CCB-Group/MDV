import React, { createContext, useContext, useState, useEffect } from 'react';

const ProjectContext = createContext(null);

export const ProjectProvider = ({ children }) => {
    // Derive initial state from window.location and URL parameters
    const { origin, pathname } = window.location;
    const flaskURL = origin + pathname;
    const urlParams = new URLSearchParams(window.location.search);
    const isPopout = urlParams.get('popout') === "true";
    const dir = urlParams.get('dir') || (isPopout ? '' : flaskURL);

    const getRoot = (dir) => {
        return dir.endsWith("/") ? dir.substring(0, dir.length - 1) : dir;
    };
    const root = getRoot(dir);
    const staticFolder = urlParams.get('static') !== null;
    const projectName = dir.split("/").pop();

    const [projectConfig, setProjectConfig] = useState({
        flaskURL,
        isPopout,
        dir,
        root,
        staticFolder,
        projectName
    });

    return (
        <ProjectContext.Provider value={{ ...projectConfig, setProjectConfig }}>
            {children}
        </ProjectContext.Provider>
    );
};

export const useProject = () => useContext(ProjectContext);
