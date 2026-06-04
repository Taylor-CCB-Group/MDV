import { useCallback, useMemo, useState } from "react";
import type { Project, ProjectAccessType, RecycledProject } from "../utils/projectUtils";
import { parseErrorResponse } from "../utils/apiUtils";
import { apiFetch, buildApiUrl } from "@/utils/mdvRouting";
import {
    projectMutationResponseSchema,
    purgeRecycleBinResponseSchema,
    recycleBinResponseSchema,
    restoreRecycledProjectResponseSchema,
} from "../utils/projectSchemas";

const useProjects = () => {
    const [projects, setProjects] = useState<Project[]>([]);
    const [isLoading, setIsLoading] = useState(false);
    const [error, setError] = useState<string | null>(null);
    const [isErrorModalOpen, setIsErrorModalOpen] = useState(false);
    const [filter, setFilter] = useState("");
    const [sortBy, setSortBy] = useState<"lastModified" | "name">(
        "lastModified",
    );
    const [sortOrder, setSortOrder] = useState<"asc" | "desc">("asc");

    const closeErrorModal = useCallback(() => {
        setIsErrorModalOpen(false);
        setError(null);
    }, []);

    const handleError = useCallback((errorMessage: string) => {
        setError(errorMessage);
        setIsErrorModalOpen(true);
    }, []);

    const fetchProjects = useCallback(async () => {
        setIsLoading(true);
        setError(null);
        try {
            const response = await apiFetch("projects", {
                headers: {
                    Accept: "application/json",
                },
            });

            if (response.ok) {
                const data: Project[] = await response.json();

                if (!Array.isArray(data)) {
                    throw new Error("Invalid response format");
                }

                const formattedProjects: Project[] = data.map((item: any) => ({
                    id: item.id === undefined || item.id === null ? "" : String(item.id),
                    name: item.name || "",
                    type: item.type === "Read-Only" ? "Read-Only" : "Editable",
                    lastModified: item.lastModified || "",
                    createdAt: item.createdAt || "",
                    owner: item.owner || [],
                    collaborators: Array.isArray(item.collaborators)
                        ? item.collaborators
                        : [],
                    numberOfStructures: item.numberOfStructures || "0",
                    numberOfImages: item.numberOfImages || "0",
                    // Assigning all permissions as true if auth is not enabled
                    permissions: item?.permissions ? {
                                read: item.permissions["can_read"],
                                edit: item.permissions["can_write"],
                                owner: item.permissions["is_owner"],
                            } : {
                                read: true,
                                edit: true,
                                owner: true,
                            },
                    thumbnail: item?.thumbnail,
                    readme: item?.readme,
                }));

                setProjects(formattedProjects);
                return;
            } else {
                const errorResponse = await parseErrorResponse({
                    response, 
                    fallbackText: "Error fetching projects. Please try again later.",
                });
                throw errorResponse;
            }
        } catch (error) {
            const errorMessage =
                error instanceof Error
                    ? error.message
                    : "Error fetching projects. Please try again later.";

            handleError(errorMessage);
            console.error("Error fetching projects:", error);
        } finally {
            setIsLoading(false);
        }
    }, [handleError]);

    const createProject = useCallback(async () => {
        setIsLoading(true);
        setError(null);

        try {
            const response = await apiFetch("create_project", {
                method: "POST",
                headers: {
                    Accept: "application/json",
                    "Content-Type": "application/json",
                },
            });

            
            if (response.ok) {
                const data = await response.json();
                // Add the new project to the local state with default values
                // These will be populated by subsequent fetch
                const newProject: Project = {
                    id: data.id,
                    name: data.name,
                    type: "Editable" as const,
                    lastModified: "",
                    createdAt: "",
                    owner: [],
                    collaborators: [],
                    numberOfStructures: "0",
                    numberOfImages: "0",
                    permissions: {
                        read: true,
                        edit: true,
                        owner: true,
                    },
                };

                setProjects((prevProjects) => [...prevProjects, newProject]);
                return newProject;
            } else {
                const errorResponse = await parseErrorResponse({
                    response, 
                    fallbackText: "Error creating project. Please try again later."
                });
                throw errorResponse;
            }

        } catch (error) {
            const errorMessage =
                error instanceof Error
                    ? error.message
                    : "Error creating project. Please try again later.";

            handleError(errorMessage);
            // todo: revisit this soon for better error handling
            throw new Error(errorMessage);
        } finally {
            setIsLoading(false);
        }
    }, [handleError]);

    const deleteProject = useCallback(
        async (id: string) => {
            setIsLoading(true);
            setError(null);

            try {
                const response = await apiFetch(`delete_project/${id}`, {
                    method: "DELETE",
                    headers: {
                        Accept: "application/json",
                    },
                });

                if (response.ok) {
                    setProjects((prevProjects) =>
                        prevProjects.filter((p) => String(p.id) !== String(id)),
                    );
                    return;
                } else {
                    const errorResponse = await parseErrorResponse({
                        response, 
                        fallbackText: "Error deleting project. Please try again later."
                    });
                    throw errorResponse;
                }
            } catch (error) {
                const errorMessage =
                    error instanceof Error
                        ? error.message
                        : "Error deleting project. Please try again later.";

                handleError(errorMessage);
                console.error("Error deleting project:", error);
            } finally {
                setIsLoading(false);
            }
        },
        [handleError],
    );

    const filteredAndSortedProjects = useMemo(() => {
        const result = projects.filter((p) =>
            p.name.toLowerCase().includes(filter.toLowerCase()),
        );

        result.sort((a, b) => {
            if (sortBy === "name") {
                return sortOrder === "asc"
                    ? a.name.localeCompare(b.name)
                    : b.name.localeCompare(a.name);
            }
            return sortOrder === "asc"
                ? new Date(a.lastModified).getTime() -
                      new Date(b.lastModified).getTime()
                : new Date(b.lastModified).getTime() -
                      new Date(a.lastModified).getTime();
        });

        return result;
    }, [projects, filter, sortBy, sortOrder]);

    const renameProject = useCallback(
        async (id: string, newName: string) => {
            setIsLoading(true);
            setError(null);

            
            try {
                if (!newName || newName.trim() === "") {
                    const errorMessage = "New name not provided";
                    handleError(errorMessage);
                    setIsLoading(false);
                    throw new Error(errorMessage);
                }
                const formData = new FormData();
                formData.append("name", newName.trim());

                const response = await apiFetch(`projects/${id}/rename`, {
                    method: "PUT",
                    body: formData,
                    headers: {
                        Accept: "application/json",
                    },
                });

                if (response.ok) {
                    setProjects((prevProjects) =>
                        prevProjects.map((project) =>
                            project.id === id
                                ? { ...project, name: newName.trim() }
                                : project,
                        ),
                    );
                    return;
                } else {
                    const errorResponse = await parseErrorResponse({
                        response, 
                        fallbackText: "Error renaming project. Please try again later."
                    });
                    throw errorResponse;
                }
            } catch (error) {
                const errorMessage =
                    error instanceof Error
                        ? error.message
                        : "Error renaming project. Please try again later.";

                handleError(errorMessage);
                console.error("Error renaming project:", error);
            } finally {
                setIsLoading(false);
            }
        },
        [handleError],
    );

    const changeProjectType = useCallback(
        async (id: string, newType: ProjectAccessType) => {
            setIsLoading(true);
            setError(null);

            try {
                const formData = new FormData();
                formData.append("type", newType.toLowerCase());

                const response = await apiFetch(`projects/${id}/access`, {
                    method: "PUT",
                    body: formData,
                    headers: {
                        Accept: "application/json",
                    },
                });

                if (response.ok) {
                    setProjects((prevProjects) =>
                        prevProjects.map((project) =>
                            project.id === id
                                ? {
                                      ...project,
                                      type: newType,
                                      lastModified: new Date().toISOString(),
                                  }
                                : project,
                        ),
                    );
                    return;
                } else {
                    const errorResponse = await parseErrorResponse({
                        response, 
                        fallbackText: "Error changing project type. Please try again later."
                    });
                    throw errorResponse;
                }
            } catch (error) {
                const errorMessage =
                    error instanceof Error
                        ? error.message
                        : "Error changing project type. Please try again later.";

                handleError(errorMessage);
                console.error("Error changing project type:", error);
            } finally {
                setIsLoading(false);
            }
        },
        [handleError],
    );

    // todo: Handle exporting larger files over 500 MB, gives error for now
    const exportProject = useCallback(
        async (id: string, name: string) => {
            setIsLoading(true);
            setError(null);

            try {
                const response = await fetch(buildApiUrl(`export_project/${id}`));
                if (response.ok) {
                    // Fetch blob from response
                    const blob = await response.blob();
                    // Create new link to download the file
                    const url = window.URL.createObjectURL(blob);
                    const link = document.createElement("a");
                    link.href = url;
                    link.download = `${name}.mdv.zip`;
                    document.body.appendChild(link);
                    link.click();
                    document.body.removeChild(link);
                    window.URL.revokeObjectURL(url);
                } else {
                    const errorResponse = await parseErrorResponse({
                        response, 
                        fallbackText: "Error exporting the project. Please try again later."
                    });
                    throw errorResponse;
                }
            } catch (error) {
                const errorMessage =
                    error instanceof Error
                        ? error.message
                        : "Error exporting the project. Please try again later.";

                handleError(errorMessage);
                console.error("Error exporting the project:", error);
            } finally {
                setIsLoading(false);
            }
    }, [handleError])

    const rescanProjects = useCallback(
        async () => {
            setIsLoading(true);
            setError(null);

            try {
                const response = await apiFetch("rescan_projects");
                if (response.ok) {
                   console.log("Rescan successful");
                   // Refresh the project list so newly discovered projects appear (works with or without auth)
                   await fetchProjects();
                } else {
                    const errorResponse = await parseErrorResponse({
                        response, 
                        fallbackText: "Error rescanning projects. Please try again later."
                    });
                    throw errorResponse;
                }
            } catch (error) {
                const errorMessage =
                    error instanceof Error
                        ? error.message
                        : "Error rescanning projects. Please try again later.";

                handleError(errorMessage);
                console.error("Error rescanning projects:", error);
            } finally {
                setIsLoading(false);
            }
        }, [handleError, fetchProjects]
    );

    const softDeleteProjects = useCallback(
        async (projectIds: string[]) => {
            setIsLoading(true);
            setError(null);
            try {
                const parsedProjectIds = projectIds.map((projectId) => Number(projectId));
                if (parsedProjectIds.some((projectId) => !Number.isInteger(projectId))) {
                    throw new Error("Invalid project id selected for deletion.");
                }

                const response = await apiFetch("projects/soft-delete", {
                    method: "POST",
                    headers: {
                        Accept: "application/json",
                        "Content-Type": "application/json",
                    },
                    body: JSON.stringify({ projectIds: parsedProjectIds }),
                });
                if (!response.ok) {
                    throw await parseErrorResponse({
                        response,
                        fallbackText: "Error moving projects to the recycle bin.",
                    });
                }
                const data = projectMutationResponseSchema.parse(await response.json());
                const deletedIds = new Set(data.deletedProjectIds.map(String));
                setProjects((previousProjects) =>
                    previousProjects.filter((project) => !deletedIds.has(String(project.id))),
                );
                return [...deletedIds];
            } catch (error) {
                const errorMessage = error instanceof Error
                    ? error.message
                    : "Error moving projects to the recycle bin.";
                handleError(errorMessage);
                throw new Error(errorMessage);
            } finally {
                setIsLoading(false);
            }
        },
        [handleError],
    );

    const fetchRecycleBin = useCallback(async () => {
        setIsLoading(true);
        setError(null);
        try {
            const response = await apiFetch("projects/recycle-bin", {
                headers: { Accept: "application/json" },
            });
            if (!response.ok) {
                throw await parseErrorResponse({
                    response,
                    fallbackText: "Error fetching the recycle bin.",
                });
            }
            return recycleBinResponseSchema.parse(await response.json());
        } catch (error) {
            const errorMessage = error instanceof Error
                ? error.message
                : "Error fetching the recycle bin.";
            handleError(errorMessage);
            return [];
        } finally {
            setIsLoading(false);
        }
    }, [handleError]);

    const purgeRecycleBin = useCallback(async () => {
        setIsLoading(true);
        setError(null);
        try {
            const response = await apiFetch("projects/recycle-bin", {
                method: "DELETE",
                headers: { Accept: "application/json" },
            });
            if (!response.ok) {
                throw await parseErrorResponse({
                    response,
                    fallbackText: "Error permanently deleting recycled projects.",
                });
            }
            const data = purgeRecycleBinResponseSchema.parse(await response.json());
            if (data.failures.length) {
                handleError(
                    `Some projects could not be permanently deleted: ${data.failures
                        .map((failure) => `${failure.projectId}: ${failure.message}`)
                        .join("; ")}`,
                );
            }
            return data;
        } catch (error) {
            const errorMessage = error instanceof Error
                ? error.message
                : "Error permanently deleting recycled projects.";
            handleError(errorMessage);
            throw new Error(errorMessage);
        } finally {
            setIsLoading(false);
        }
    }, [handleError]);

    const restoreRecycledProject = useCallback(
        async (projectId: string) => {
            setIsLoading(true);
            setError(null);
            try {
                const response = await apiFetch(`projects/recycle-bin/${projectId}/restore`, {
                    method: "POST",
                    headers: { Accept: "application/json" },
                });
                if (!response.ok) {
                    throw await parseErrorResponse({
                        response,
                        fallbackText: "Error restoring recycled project.",
                    });
                }
                await fetchProjects();
                return restoreRecycledProjectResponseSchema.parse(await response.json());
            } catch (error) {
                const errorMessage = error instanceof Error
                    ? error.message
                    : "Error restoring recycled project.";
                handleError(errorMessage);
                throw new Error(errorMessage);
            } finally {
                setIsLoading(false);
            }
        },
        [fetchProjects, handleError],
    );

    return {
        projects: filteredAndSortedProjects,
        isLoading,
        error,
        isErrorModalOpen,
        closeErrorModal,
        fetchProjects,
        createProject,
        deleteProject,
        softDeleteProjects,
        fetchRecycleBin,
        purgeRecycleBin,
        restoreRecycledProject,
        renameProject,
        changeProjectType,
        setFilter,
        setSortBy,
        setSortOrder,
        sortBy,
        exportProject,
        rescanProjects,
    };
};

export default useProjects;
