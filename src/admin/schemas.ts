import { z } from "zod/v4";

export const adminPermissionSchema = z.enum(["view", "edit", "owner"]);

export const adminUserSchema = z.object({
    id: z.number(),
    email: z.string(),
    firstName: z.string(),
    lastName: z.string(),
    isActive: z.boolean(),
    isAdmin: z.boolean(),
});

export const adminProjectSchema = z.object({
    id: z.number(),
    name: z.string(),
    path: z.string(),
    accessLevel: z.string(),
    isPublic: z.boolean(),
    isDeleted: z.boolean(),
    updatedAt: z.string().nullable(),
});

export const adminSessionSchema = z.object({
    authEnabled: z.boolean(),
    isAdmin: z.boolean(),
    permissions: z.array(z.string()),
    user: z.object({
        id: z.number(),
        email: z.string(),
        is_admin: z.boolean().optional(),
        synthetic: z.boolean().optional(),
    }),
});

export const adminProjectAccessSchema = z.object({
    projectId: z.number(),
    permission: adminPermissionSchema,
    canRead: z.boolean(),
    canWrite: z.boolean(),
    isOwner: z.boolean(),
});

export const createAdminUserPayloadSchema = z.object({
    email: z.string(),
    firstName: z.string().default(""),
    lastName: z.string().default(""),
    projectAccess: z.array(z.object({
        projectId: z.number(),
        permission: adminPermissionSchema,
    })).default([]),
});

export const createAdminUserResultSchema = z.object({
    user: adminUserSchema,
    projectAccess: z.array(adminProjectAccessSchema),
    created: z.boolean(),
});

export const adminProjectMemberSchema = z.object({
    user: adminUserSchema,
    projectAccess: adminProjectAccessSchema,
});

export const addProjectMemberPayloadSchema = z.object({
    userId: z.number(),
    permission: adminPermissionSchema,
});

export const updateProjectMemberPayloadSchema = z.object({
    permission: adminPermissionSchema,
});

export const adminUsersResponseSchema = z.object({
    users: z.array(adminUserSchema),
});

export const adminProjectsResponseSchema = z.object({
    projects: z.array(adminProjectSchema),
});

export const adminProjectMembersResponseSchema = z.object({
    members: z.array(adminProjectMemberSchema),
});

export const adminErrorResponseSchema = z.object({
    error: z.string(),
});

export type AdminPermission = z.infer<typeof adminPermissionSchema>;
export type AdminUser = z.infer<typeof adminUserSchema>;
export type AdminProject = z.infer<typeof adminProjectSchema>;
export type AdminSession = z.infer<typeof adminSessionSchema>;
export type CreateAdminUserPayload = z.infer<typeof createAdminUserPayloadSchema>;
export type CreateAdminUserResult = z.infer<typeof createAdminUserResultSchema>;
export type AdminProjectMember = z.infer<typeof adminProjectMemberSchema>;
