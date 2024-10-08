import { useState, useCallback, useMemo } from 'react';

interface Project {
    id: string;
    name: string;
    type: "Editable" | "Read-Only";
    lastModified: string;
    createdAt: string;
    owner: string;
    collaborators: string[];
    numberOfStructures: string;
    numberOfImages: string;
}

const useProjects = () => {
  const [projects, setProjects] = useState<Project[]>([]);
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [filter, setFilter] = useState("");
  const [sortBy, setSortBy] = useState<"lastModified" | "name">("lastModified");
  const [sortOrder, setSortOrder] = useState<"asc" | "desc">("desc");

  const fetchProjects = useCallback(async () => {
    setIsLoading(true);
    setError(null);
    try {
      const response = await fetch("projects");
      if (!response.ok) {
        throw new Error('Failed to fetch projects');
      }
      const data: Project[] = await response.json();
      const formattedProjects: Project[] = data.map((item: any) => ({
        id: item.id || '',
        name: item.name || '',
        type: item.type || '',
        lastModified: item.lastModified || '',
        createdAt: item.createdAt || '',
        owner: item.owner || '',
        collaborators: item.collaborators || [],
        numberOfStructures: item.numberOfStructures || '0',
        numberOfImages: item.numberOfImages || '0',
      }));
      setProjects(formattedProjects);
    } catch (error) {
      setError('Error fetching projects. Please try again later.');
      console.error("Error fetching projects:", error);
    } finally {
      setIsLoading(false);
    }
  }, []);

  const createProject = useCallback(async (projectName: string) => {
    setIsLoading(true);
    setError(null);
    try {
      const response = await fetch("create_project", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
        },
        body: JSON.stringify({
          id: projectName, // will be ignored with real backend, but currently used in mdv_desktop
        }),
      });

      if (!response.ok) {
        throw new Error(response.statusText);
      }

      const data = await response.json();
      setProjects(prevProjects => [...prevProjects, data]);
      return data;
      
    } catch (error) {
      setError('Error creating project. Please try again later.');
      console.error("Error creating project:", error);
    } finally {
      setIsLoading(false);
    }
  }, []);

  const deleteProject = useCallback(async (id: string) => {
    try {
      const response = await fetch(`delete_project/${id}`, {
        method: "DELETE",
      });
      if (!response.ok) {
        throw new Error("Failed to delete project");
      }
      setProjects((prevProjects) => prevProjects.filter((p) => p.id !== id));
    } catch (error) {
      setError('Error deleting project. Please try again later.');
      console.error("Error deleting project:", error);
    }
  }, []);

  const filteredAndSortedProjects = useMemo(() => {
    const result = projects.filter((p) =>
      p.name.toLowerCase().includes(filter.toLowerCase())
    );

    result.sort((a, b) => {
      if (sortBy === "name") {
        return sortOrder === "asc"
          ? a.name.localeCompare(b.name)
          : b.name.localeCompare(a.name);
      }
      return sortOrder === "asc"
        ? new Date(a.lastModified).getTime() - new Date(b.lastModified).getTime()
        : new Date(b.lastModified).getTime() - new Date(a.lastModified).getTime();
    });

    return result;
  }, [projects, filter, sortBy, sortOrder]);

  return {
    projects: filteredAndSortedProjects,
    isLoading,
    error,
    fetchProjects,
    createProject,
    deleteProject,
    setFilter,
    setSortBy,
    setSortOrder,
    sortBy,
  };
};

export default useProjects;