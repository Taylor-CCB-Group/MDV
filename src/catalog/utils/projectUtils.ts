export type SortBy = "lastModified" | "name";
export type SortOrder = "asc" | "desc";
// todo - consider relating this to what is in state.json (used in ChartManager)
export type ProjectAccessType = "Editable" | "Read-Only";
export type Permissions = {
    read: boolean;
    edit: boolean;
    owner: boolean;
};
export interface Project {
    id: string;
    name: string;
    type: ProjectAccessType;
    lastModified: string;
    createdAt: string;
    owner: string[];
    collaborators: string[];
    numberOfStructures: string;
    numberOfImages: string;
    permissions: Permissions;
    thumbnail?: string;
    readme?: string;
}

export const sortProjects = (
    projects: Project[],
    sortBy: SortBy,
    sortOrder: SortOrder = sortBy === "lastModified" ? "desc" : "desc"
): Project[] => {
    return [...projects].sort((a, b) => {
        let comparison = 0;
        
        if (sortBy === "name") {
            comparison = b.name.localeCompare(a.name, undefined, { sensitivity: 'base' }); // Descending by default
        } else if (sortBy === "lastModified") {
            comparison = new Date(b.lastModified).getTime() - new Date(a.lastModified).getTime(); // Most recent first
        }
        
        return sortOrder === "asc" ? -comparison : comparison;
    });
};
