import { useCallback, useMemo, useState } from "react";
import type { Project, ProjectAccessType } from "../utils/projectUtils";
import { parseErrorResponse } from "../utils/apiUtils";

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
            const response = await fetch("projects", {
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
                    id: item.id || "",
                    name: item.name || "",
                    type: item.type === "Read-Only" ? "Read-Only" : "Editable",
                    lastModified: item.lastModified || "",
                    createdAt: item.createdAt || "",
                    owner: item.owner || "",
                    collaborators: Array.isArray(item.collaborators)
                        ? item.collaborators
                        : [],
                    numberOfStructures: item.numberOfStructures || "0",
                    numberOfImages: item.numberOfImages || "0",
                    thumbnail: item?.thumbnail,
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
            const response = await fetch("create_project", {
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
                    owner: "",
                    collaborators: [],
                    numberOfStructures: "0",
                    numberOfImages: "0",
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
        } finally {
            setIsLoading(false);
        }
    }, [handleError]);

    const deleteProject = useCallback(
        async (id: string) => {
            setIsLoading(true);
            setError(null);

            try {
                const response = await fetch(`delete_project/${id}`, {
                    method: "DELETE",
                    headers: {
                        Accept: "application/json",
                    },
                });

                if (response.ok) {
                    setProjects((prevProjects) =>
                        prevProjects.filter((p) => p.id !== id),
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

                const response = await fetch(`projects/${id}/rename`, {
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

                const response = await fetch(`projects/${id}/access`, {
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

    return {
        projects: filteredAndSortedProjects,
        isLoading,
        error,
        isErrorModalOpen,
        closeErrorModal,
        fetchProjects,
        createProject,
        deleteProject,
        renameProject,
        changeProjectType,
        setFilter,
        setSortBy,
        setSortOrder,
        sortBy,
    };
};

export default useProjects;
