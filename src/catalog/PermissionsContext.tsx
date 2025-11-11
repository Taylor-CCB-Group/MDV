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

// Default state
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

interface PermissionsContextType {
    permissions: ProjectOperationPermissions;
    isLoading: boolean;
    error: string | null;
    isProjectManagerExists: boolean;
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
    const [isProjectManagerExists, setIsProjectManagerExists] = useState(false);

    const fetchPermissions = useCallback(async () => {
        setIsLoading(true);
        setError(null);
        try {
            const response = await fetch("extension_config");
            if (!response.ok) {
                // Fallback back to default permissions
                console.warn(
                    `Permissions endpoint not found (${response.status}). Falling back to default permissions.`,
                );
                setPermissions(defaultPermissions);
                return;
            }
            const data = await response.json();
            if (data.project_manager) {
                setIsProjectManagerExists(true);
            }

            const projectPermissions: ProjectOperationPermissions = {
                createProject: data.project_manager.createProject ?? true,
                importProject: data.project_manager.importProject ?? true,
                deleteProject: data.project_manager.deleteProject ?? true,
                renameProject: data.project_manager.renameProject ?? true,
                changeProjectAccess: data.project_manager.changeProjectAccess ?? true,
                exportProject: data.project_manager.exportProject ?? true,
                shareProject: data.project_manager.shareProject ?? true,
                editUserPermissions: data.project_manager.editUserPermissions ?? true,
                removeUserFromProject: data.project_manager.removeUserFromProject ?? true,
            };
            setPermissions(projectPermissions);
        } catch (e) {
            console.error("Failed to fetch project permissions:", e);
            setError("Failed to load permissions.");
            // Fallback to default permissions on error
            setPermissions(defaultPermissions);
            setIsProjectManagerExists(false);
        } finally {
            setIsLoading(false);
        }
    }, []);

    useEffect(() => {
        fetchPermissions();
    }, [fetchPermissions]);

    const value = { permissions, isLoading, error, isProjectManagerExists };

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
