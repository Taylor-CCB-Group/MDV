import {
    type ReactNode,
    createContext,
    useCallback,
    useContext,
    useEffect,
    useState,
} from "react";

export interface ProjectOperationPermissions {
    createProject: boolean;
    importProject: boolean;
    deleteProject: boolean;
    renameProject: boolean;
    changeProjectAccess: boolean;
    exportProject: boolean;
    shareProject: boolean;
    editUserPermissions: boolean;
    removeUserFromProject: boolean;
}

// Default state (All false)
const defaultPermissions: ProjectOperationPermissions = {
    createProject: false,
    importProject: false,
    deleteProject: false,
    renameProject: false,
    changeProjectAccess: false,
    exportProject: false,
    shareProject: false,
    editUserPermissions: false,
    removeUserFromProject: false,
};

const enabledPermissions: ProjectOperationPermissions = {
    createProject: true,
    importProject: true,
    deleteProject: true,
    renameProject: true,
    changeProjectAccess: true,
    exportProject: true,
    shareProject: true,
    editUserPermissions: true,
    removeUserFromProject: true,
};

interface PermissionsContextType {
    permissions: ProjectOperationPermissions;
    isLoading: boolean;
    error: string | null;
    isPublicPage: boolean;
}

const PermissionsContext = createContext<PermissionsContextType | undefined>(
    undefined,
);

/**
 * Provider component that fetches project operation permissions
 */
export const PermissionsProvider: React.FC<{ children: ReactNode }> = ({
    children,
}) => {
    const [permissions, setPermissions] =
        useState<ProjectOperationPermissions>(defaultPermissions);
    const [isLoading, setIsLoading] = useState(true);
    const [error, setError] = useState<string | null>(null);
    //! Public page is true by default
    const [isPublicPage, setIsPublicPage] = useState(true);

    const fetchPermissions = useCallback(async () => {
        setIsLoading(true);
        setError(null);
        try {
            const response = await fetch("extension_config");
            if (!response.ok) {
                if (response.status === 404) {
                    setPermissions(enabledPermissions);
                    setIsPublicPage(false);
                    return;
                }
                // Fallback back to default permissions
                console.warn(
                    `Permissions endpoint not found (${response.status}). Falling back to default permissions.`,
                );
                setPermissions(defaultPermissions);
                return;
            }
            const data = await response.json();

            if (data.project_manager) {
                const project_manager = data.project_manager;

                const pmOperations = Object.values(project_manager).every((operation) => !operation);
                // Public page is false if any one of the values is true, else it's true by default
                if (!pmOperations) {
                    setIsPublicPage(false);
                }
                
                // If any permission doesn't exist, then we assign it as false
                const projectPermissions: ProjectOperationPermissions = {
                    createProject: data.project_manager?.createProject ?? false,
                    importProject: data.project_manager?.importProject ?? false,
                    deleteProject: data.project_manager?.deleteProject ?? false,
                    renameProject: data.project_manager?.renameProject ?? false,
                    changeProjectAccess: data.project_manager?.changeProjectAccess ?? false,
                    exportProject: data.project_manager?.exportProject ?? false,
                    shareProject: data.project_manager?.shareProject ?? false,
                    editUserPermissions: data.project_manager?.editUserPermissions ?? false,
                    removeUserFromProject: data.project_manager?.removeUserFromProject ?? false,
                };
                setPermissions(projectPermissions);
            } else {
                setPermissions(enabledPermissions);
                setIsPublicPage(false);
            }
        } catch (e) {
            console.error("Failed to fetch project permissions:", e);
            setError("Failed to load permissions.");
            setPermissions(enabledPermissions);
            setIsPublicPage(false);
        } finally {
            setIsLoading(false);
        }
    }, []);

    useEffect(() => {
        fetchPermissions();
    }, [fetchPermissions]);

    const value = { permissions, isLoading, error, isPublicPage };

    return (
        <PermissionsContext.Provider value={value}>
            {children}
        </PermissionsContext.Provider>
    );
};

// A custom hook to consume Project Permissions
const usePermissions = (): PermissionsContextType => {
    const context = useContext(PermissionsContext);
    if (context === undefined) {
        throw new Error("usePermissions must be used within a PermissionsProvider");
    }
    return context;
};

export default usePermissions;
