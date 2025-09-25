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
    changeProjectType: boolean;
    exportProject: boolean;
    shareProject: boolean;
}

// Default state
const defaultPermissions: ProjectOperationPermissions = {
    createProject: false,
    importProject: false,
    deleteProject: false,
    renameProject: false,
    changeProjectType: false,
    exportProject: false,
    shareProject: false,
};

interface PermissionsContextType {
    permissions: ProjectOperationPermissions;
    isLoading: boolean;
    error: string | null;
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
            const projectPermissions: ProjectOperationPermissions = {
                createProject: data.projectManager.createProject ?? true,
                importProject: data.projectManager.importProject ?? true,
                deleteProject: data.projectManager.deleteProject ?? true,
                renameProject: data.projectManager.renameProject ?? true,
                changeProjectType: data.projectManager.changeProjectType ?? true,
                exportProject: data.projectManager.exportProject ?? true,
                shareProject: data.projectManager.shareProject ?? true,
            };
            setPermissions(projectPermissions);
        } catch (e) {
            console.error("Failed to fetch project permissions:", e);
            setError("Failed to load permissions.");
            // Fallback to default permissions on error
            setPermissions(defaultPermissions);
        } finally {
            setIsLoading(false);
        }
    }, []);

    useEffect(() => {
        fetchPermissions();
    }, [fetchPermissions]);

    const value = { permissions, isLoading, error };

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
