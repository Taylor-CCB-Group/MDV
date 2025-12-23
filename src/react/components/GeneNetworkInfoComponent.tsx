import { useEffect, useRef, useState } from "react";
import { useQuery } from "@tanstack/react-query";
import { Paper, Typography, Link } from "@mui/material";
import z from "zod";

// Zod schema for GeneNetwork API response
const geneNetworkApiResponseSchema = z.object({
    comment: z.string().optional(),
    version: z.string().optional(),
    gene: z.object({
        id: z.string(),
        index_: z.number().optional(),
        name: z.string(),
        biotype: z.string().optional(),
        chr: z.string().optional(),
        start: z.number().optional(),
        stop: z.number().optional(),
        strand: z.number().optional(),
        description: z.string(),
        biomartRelease: z.string().optional(),
        assemblyRelease: z.string().optional(),
        transcripts: z.array(z.string()).optional(),
        genePredScore: z.string().optional(),
    }),
    databases: z.array(z.object({
        id: z.string(),
        name: z.string(),
        fullName: z.string(),
        description: z.string(),
        url: z.string(),
        inMemory: z.boolean(),
    })).optional(),
    pathways: z.object({
        annotated: z.array(z.object({
            href: z.string(),
            zScore: z.number(),
            pValue: z.number(),
        })).optional(),
        predicted: z.array(z.object({
            href: z.string(),
            zScore: z.number(),
            pValue: z.number(),
        })).optional(),
    }).optional(),
    celltypes: z.object({
        fixed: z.object({
            header: z.array(z.any()).optional(),
            indices: z.record(z.string(), z.number()).optional(),
        }).optional(),
        values: z.object({
            avg: z.array(z.number()).optional(),
            stdev: z.array(z.number()).optional(),
            z: z.array(z.number()).optional(),
            auc: z.array(z.number()).optional(),
        }).optional(),
        transcriptBars: z.record(z.string(), z.array(z.number().nullable())).optional(),
    }).optional(),
}).passthrough();

export type GeneNetworkApiResponse = z.infer<typeof geneNetworkApiResponseSchema>;

// Mapped interface for component use
export interface GeneInfo {
    geneId: string;
    geneName: string;
    geneDescription: string;
    geneProperties?: {
        alias?: string[];
    };
}

interface GeneNetworkInfoComponentProps {
    geneId: string;
}

/**
 * Component that fetches and displays gene information from GeneNetwork.nl.
 * Uses lazy loading with IntersectionObserver to only query when visible.
 */
export function GeneNetworkInfoComponent({ geneId }: GeneNetworkInfoComponentProps) {
    const containerRef = useRef<HTMLDivElement>(null);
    const [isVisible, setIsVisible] = useState(false);

    // Set up IntersectionObserver for lazy loading
    useEffect(() => {
        if (!containerRef.current) return;

        const observer = new IntersectionObserver(
            (entries) => {
                if (entries[0].isIntersecting && !isVisible) {
                    setIsVisible(true);
                }
            },
            { rootMargin: "100px" }, // Start loading slightly before it's visible
        );

        observer.observe(containerRef.current);

        return () => {
            observer.disconnect();
        };
    }, [isVisible]);

    // Only query when visible and geneId is provided
    const { data: geneInfo, isLoading, error } = useQuery<GeneInfo>({
        queryKey: ["geneInfo", geneId],
        queryFn: async () => {
            if (!geneId) throw new Error("No gene ID provided");
            const response = await fetch(
                `https://www.genenetwork.nl/api/v1/gene/${geneId}`,
            );
            if (!response.ok) {
                const errorData = await response.json().catch(() => null);
                throw new Error(
                    errorData?.message ||
                        `HTTP error! status: ${response.status}`,
                );
            }
            const jsonData = await response.json();
            // Parse and validate with Zod
            const parsed = geneNetworkApiResponseSchema.parse(jsonData);
            // Map API response to component's expected format
            return {
                geneId: parsed.gene.id,
                geneName: parsed.gene.name,
                geneDescription: parsed.gene.description,
                geneProperties: {
                    // Aliases not present in API response, but keeping structure for future use
                    alias: undefined,
                },
            };
        },
        enabled: isVisible && !!geneId,
        // staleTime: 5 * 60 * 1000, // Cache for 5 minutes
    });

    return (
        <Paper
            ref={containerRef}
            variant="outlined"
            sx={{
                my: 1,
                p: 1.5,
                borderRadius: 1,
                minHeight: 100,
            }}
        >
            {!isVisible && (
                <Typography variant="body2" color="text.secondary">
                    Scroll to load gene info...
                </Typography>
            )}
            {isVisible && isLoading && (
                <Typography variant="body2">
                    Loading gene info for {geneId}...
                </Typography>
            )}
            {isVisible && error && (
                <Typography variant="body2" color="error">
                    Error fetching gene info:{" "}
                    {error instanceof Error ? error.message : "Unknown error"}
                </Typography>
            )}
            {isVisible && geneInfo && (
                <>
                    <Typography variant="subtitle1" fontWeight={600} gutterBottom>
                        Gene Information for {geneInfo.geneName}
                    </Typography>
                    <Typography variant="body2">
                        <strong>Gene ID:</strong> {geneInfo.geneId}
                    </Typography>
                    <Typography variant="body2" sx={{ mt: 0.5 }}>
                        <strong>Description:</strong> {geneInfo.geneDescription}
                    </Typography>
                    {geneInfo.geneProperties?.alias &&
                        geneInfo.geneProperties.alias.length > 0 && (
                            <Typography variant="body2" sx={{ mt: 0.5 }}>
                                <strong>Aliases:</strong>{" "}
                                {geneInfo.geneProperties.alias.join(", ")}
                            </Typography>
                        )}
                    <Typography variant="body2" sx={{ mt: 1 }}>
                        <Link
                            href={`https://www.genenetwork.nl/gene/${geneInfo.geneId}`}
                            target="_blank"
                            rel="noopener noreferrer"
                        >
                            More info on GeneNetwork
                        </Link>
                    </Typography>
                </>
            )}
        </Paper>
    );
}

