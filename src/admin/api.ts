export type AdminSession = {
    authEnabled: boolean;
    isAdmin: boolean;
    permissions: string[];
    user: {
        id: number;
        email: string;
        is_admin?: boolean;
        synthetic?: boolean;
    };
};

export type AdminUser = {
    id: number;
    email: string;
    firstName: string;
    lastName: string;
    isActive: boolean;
    isAdmin: boolean;
};

export type AdminProject = {
    id: number;
    name: string;
    path: string;
    accessLevel: string;
    isPublic: boolean;
    isDeleted: boolean;
    updatedAt: string | null;
};

async function requestJson<T>(path: string, init?: RequestInit): Promise<T> {
    const response = await fetch(path, {
        credentials: "same-origin",
        ...init,
    });

    if (!response.ok) {
        const body = await response.json().catch(() => undefined);
        const message =
            body && typeof body.error === "string"
                ? body.error
                : `Request failed with ${response.status}`;
        throw new Error(message);
    }

    return response.json();
}

export const adminApi = {
    session: () => requestJson<AdminSession>("/admin/api/session"),
    users: () => requestJson<{ users: AdminUser[] }>("/admin/api/users"),
    projects: () =>
        requestJson<{ projects: AdminProject[] }>("/admin/api/projects"),
};
