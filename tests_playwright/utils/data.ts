//! This needs to be in sync with the return type of the /projects api
// maybe make use of zod schema in future
export type MockProject = {
    id: string;
    name: string;
    type: "Editable" | "Read-Only";
    lastModified: string;
    createdAt?: string;
    owner?: any[];
    collaborators?: any[];
    numberOfStructures?: string;
    numberOfImages?: string;
    permissions?: { can_read: boolean; can_write: boolean; is_owner: boolean };
    readme?: string;
};

export const project = (overrides: Partial<MockProject> = {}): MockProject => ({
    id: "p1",
    name: "Alpha",
    type: "Editable",
    lastModified: "2025-01-01T00:00:00Z",
    createdAt: "2025-01-01T00:00:00Z",
    owner: [],
    collaborators: [],
    numberOfStructures: "0",
    numberOfImages: "0",
    permissions: { can_read: true, can_write: true, is_owner: true },
    readme: "Readme",
    ...overrides,
});
