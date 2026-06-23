import type { z } from "zod/v4";
import {
    addProjectMemberPayloadSchema,
    adminErrorResponseSchema,
    adminProjectMemberSchema,
    adminProjectMembersResponseSchema,
    adminProjectsResponseSchema,
    adminSessionSchema,
    adminUsersResponseSchema,
    createAdminUserPayloadSchema,
    createAdminUserResultSchema,
    updateProjectMemberPayloadSchema,
    type AdminPermission,
    type AdminProject,
    type AdminProjectMember,
    type AdminSession,
    type AdminUser,
    type CreateAdminUserPayload,
    type CreateAdminUserResult,
} from "./schemas";

export type {
    AdminPermission,
    AdminProject,
    AdminProjectMember,
    AdminSession,
    AdminUser,
    CreateAdminUserPayload,
    CreateAdminUserResult,
};

async function requestJson<T>(
    path: string,
    schema: z.ZodType<T>,
    init?: RequestInit,
): Promise<T> {
    const response = await fetch(path, {
        credentials: "same-origin",
        ...init,
    });

    if (!response.ok) {
        const body = await response.json().catch(() => undefined);
        const parsed = adminErrorResponseSchema.safeParse(body);
        const message = parsed.success
            ? parsed.data.error
            : `Request failed with ${response.status}`;
        throw new Error(message);
    }

    return schema.parse(await response.json());
}

export const adminApi = {
    session: () => requestJson<AdminSession>("/admin/api/session", adminSessionSchema),
    users: () => requestJson("/admin/api/users", adminUsersResponseSchema),
    createUser: (body: CreateAdminUserPayload) => {
        const payload = createAdminUserPayloadSchema.parse(body);
        return requestJson<CreateAdminUserResult>("/admin/api/users", createAdminUserResultSchema, {
            method: "POST",
            headers: {
                "Content-Type": "application/json",
            },
            body: JSON.stringify(payload),
        });
    },
    projects: () =>
        requestJson("/admin/api/projects", adminProjectsResponseSchema),
    projectMembers: (projectId: number) =>
        requestJson(
            `/admin/api/projects/${projectId}/users`,
            adminProjectMembersResponseSchema,
        ),
    addProjectMember: (
        projectId: number,
        body: { userId: number; permission: AdminPermission },
    ) => {
        const payload = addProjectMemberPayloadSchema.parse(body);
        return requestJson<AdminProjectMember>(
            `/admin/api/projects/${projectId}/users`,
            adminProjectMemberSchema,
            {
                method: "POST",
                headers: {
                    "Content-Type": "application/json",
                },
                body: JSON.stringify(payload),
            },
        );
    },
    updateProjectMember: (
        projectId: number,
        userId: number,
        body: { permission: AdminPermission },
    ) => {
        const payload = updateProjectMemberPayloadSchema.parse(body);
        return requestJson<AdminProjectMember>(
            `/admin/api/projects/${projectId}/users/${userId}`,
            adminProjectMemberSchema,
            {
                method: "PATCH",
                headers: {
                    "Content-Type": "application/json",
                },
                body: JSON.stringify(payload),
            },
        );
    },
    removeProjectMember: (projectId: number, userId: number) =>
        fetch(`/admin/api/projects/${projectId}/users/${userId}`, {
            method: "DELETE",
            credentials: "same-origin",
        }).then((response) => {
            if (!response.ok) {
                return response.json().then((body: unknown) => {
                    const parsed = adminErrorResponseSchema.safeParse(body);
                    const message = parsed.success
                        ? parsed.data.error
                        : `Request failed with ${response.status}`;
                    throw new Error(message);
                });
            }
        }),
};
