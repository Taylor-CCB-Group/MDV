import { useEffect, useRef, useState, type MouseEvent, type KeyboardEvent } from "react";
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
}

interface GeneNetworkInfoComponentProps {
    geneId: string;
    highlightCount?: number;
    isHighlighted?: boolean;
    onCardClick?: (event: MouseEvent<HTMLDivElement>) => void;
    onCardKeyDown?: (event: KeyboardEvent<HTMLDivElement>) => void;
}

/**
 * Component that fetches and displays gene information from GeneNetwork.nl.
 * Uses lazy loading with IntersectionObserver to only query when visible.
 */
export function GeneNetworkInfoComponent({
    geneId,
    isHighlighted = false,
    onCardClick,
    onCardKeyDown,
}: GeneNetworkInfoComponentProps) {
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
            };
        },
        enabled: isVisible && !!geneId,
        retry: false,
        // staleTime: 5 * 60 * 1000, // Cache for 5 minutes
    });

    const handleCardClick = (event: MouseEvent<HTMLDivElement>) => {
        // Don't trigger selection when clicking on links within the card
        const target = event.target as HTMLElement | null;
        if (target?.closest("a")) {
            return;
        }
        onCardClick?.(event);
    };

    const handleKeyDown = (event: KeyboardEvent<HTMLDivElement>) => {
        // Handle Space key for selection (prevent default scrolling)
        if (event.key === " " || event.key === "Spacebar") {
            event.preventDefault();
            onCardKeyDown?.(event);
            return;
        }
        // Allow other keys to bubble up (e.g., Tab for navigation)
        onCardKeyDown?.(event);
    };

    return (
        <Paper
            ref={containerRef}
            variant={isHighlighted ? "elevation" : "outlined"}
            elevation={isHighlighted ? 2 : 0}
            tabIndex={onCardClick || onCardKeyDown ? 0 : undefined}
            sx={{
                p: 1.5,
                m: 1,
                borderRadius: 1,
                minHeight: 100,
                border: isHighlighted ? "2px solid" : undefined,
                borderColor: isHighlighted ? "primary.main" : undefined,
                cursor: onCardClick ? "pointer" : "default",
                outline: "none",
                "&:focus-visible": {
                    outline: "2px solid",
                    outlineColor: "primary.main",
                    outlineOffset: "2px",
                },
            }}
            onClick={handleCardClick}
            onKeyDown={handleKeyDown}
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
                <>
                <Typography variant="body2" color="error">
                    Error fetching gene info:{" "}
                    {error instanceof Error ? error.message : "Unknown error"}
                </Typography>
                <Link
                    href={`https://www.genecards.org/cgi-bin/carddisp.pl?gene=${geneId}`}
                    target="_blank"
                    rel="noopener noreferer"
                    tabIndex={-1}
                >
                    Try opening on genecards.org
                </Link>
                </>
            )}
            {isVisible && geneInfo && (
                <>
                    <div className="flex items-start justify-between">
                        <Typography variant="subtitle1" fontWeight={600} gutterBottom>
                            Gene Information for {geneInfo.geneName}
                        </Typography>
                    </div>
                    <Typography variant="body2">
                        <strong>Gene ID:</strong> {geneInfo.geneId}
                    </Typography>
                    <Typography variant="body2" sx={{ mt: 0.5 }}>
                        <strong>Description:</strong> {geneInfo.geneDescription}
                    </Typography>
                    <Typography variant="body2" sx={{ mt: 1 }}>
                        <Link
                            href={`https://www.genenetwork.nl/gene/${geneInfo.geneId}`}
                            target="_blank"
                            rel="noopener noreferrer"
                            tabIndex={-1}
                        >
                            More info on GeneNetwork
                        </Link>
                        {" - "}
                        <Link
                            href={`https://www.genecards.org/cgi-bin/carddisp.pl?gene=${geneId}`}
                            target="_blank"
                            rel="noopener noreferer"
                            tabIndex={-1}
                        >
                            Open on genecards.org
                        </Link>
                    </Typography>
                </>
            )}
        </Paper>
    );
}

