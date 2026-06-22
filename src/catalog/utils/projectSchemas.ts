import { z } from "zod";

export const projectMutationResponseSchema = z.object({
    deletedProjectIds: z.array(z.number()),
});

export const recycledProjectSchema = z.object({
    id: z.union([z.number(), z.string()]).transform(String),
    name: z.string(),
    deletedTimestamp: z.string().nullable().optional().transform((value) => value ?? undefined),
});

export const recycleBinResponseSchema = z.array(recycledProjectSchema);

export const purgeRecycleBinResponseSchema = z.object({
    deletedProjectIds: z.array(z.number()),
    failures: z.array(z.object({
        projectId: z.number(),
        message: z.string(),
    })),
});

export const restoreRecycledProjectResponseSchema = z.object({
    restoredProjectId: z.number(),
});

export type RecycledProject = z.infer<typeof recycledProjectSchema>;
export type PurgeRecycleBinResponse = z.infer<typeof purgeRecycleBinResponseSchema>;
export type RestoreRecycledProjectResponse = z.infer<typeof restoreRecycledProjectResponseSchema>;
