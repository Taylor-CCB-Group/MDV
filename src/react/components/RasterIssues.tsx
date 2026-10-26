import type { SpatialData } from "@spatialdata/core";
import type { RenderStack } from "@spatialdata/layers";
import { useEffect, useState } from "react";
import {
    assessRasterLevels,
    describeRasterIssue,
    loadRasterLevels,
    RASTER_REENCODE_ADVICE,
    type RasterIssue,
} from "../spatialdata/raster_health";

type RasterElementType = "image" | "labels";
export type ElementIssues = { elementType: RasterElementType; elementKey: string; issues: RasterIssue[] };
const NO_ISSUES: ElementIssues[] = [];

/**
 * `"image:key|labels:key"` for the raster entries in the stack. Entries are mutated in place
 * for cosmetic edits, so this is the value the metadata fetch is keyed on — not the array.
 * Call it from an observer so stack edits are tracked.
 */
export function rasterEntriesKey(stack: RenderStack | undefined) {
    const keys = new Set<string>();
    for (const entry of stack?.entries ?? []) {
        if (entry.kind !== "spatial") continue;
        const { elementType, elementKey } = entry.source;
        if (elementType === "image" || elementType === "labels") keys.add(`${elementType}:${elementKey}`);
    }
    return [...keys].sort().join("|");
}

/**
 * Image and labels elements in the render stack that are laid out in a way that makes them
 * slow to load and render. See `raster_health.ts`.
 */
export function useRasterIssues(spatialData: SpatialData | null | undefined, entriesKey: string) {
    const [result, setResult] = useState<{ spatialData: SpatialData; entriesKey: string; issues: ElementIssues[] }>();
    useEffect(() => {
        if (!spatialData || !entriesKey) return;
        let cancelled = false;
        Promise.all(
            entriesKey.split("|").map(async (key): Promise<ElementIssues | null> => {
                const separator = key.indexOf(":");
                const elementType = key.slice(0, separator);
                const elementKey = key.slice(separator + 1);
                if (elementType !== "image" && elementType !== "labels") return null;
                const element =
                    elementType === "image" ? spatialData.images?.[elementKey] : spatialData.labels?.[elementKey];
                if (!element) return null;
                try {
                    const found = assessRasterLevels(await loadRasterLevels(element));
                    return found.length ? { elementType, elementKey, issues: found } : null;
                } catch (error) {
                    console.warn(`could not read array metadata for ${elementType} '${elementKey}'`, error);
                    return null;
                }
            }),
        ).then((results) => {
            if (cancelled) return;
            const issues = results.filter((item): item is ElementIssues => item !== null);
            setResult({ spatialData, entriesKey, issues });
        });
        return () => {
            cancelled = true;
        };
    }, [spatialData, entriesKey]);
    // A result for a different store or stack is stale until the new fetch lands.
    return result && result.spatialData === spatialData && result.entriesKey === entriesKey ? result.issues : NO_ISSUES;
}

/** Which raster elements are slow to load, and why, with the advice for fixing them. */
export function RasterIssueList({ issues }: { issues: ElementIssues[] }) {
    return (
        <div className="space-y-1">
            {issues.map((item) => (
                <div key={`${item.elementType}:${item.elementKey}`}>
                    <div className="font-medium">
                        {item.elementType === "labels" ? "Labels" : "Image"} '{item.elementKey}'
                    </div>
                    <ul className="ml-4 list-disc">
                        {item.issues.map((issue) => (
                            <li key={issue.kind}>{describeRasterIssue(issue)}</li>
                        ))}
                    </ul>
                </div>
            ))}
            <p className="opacity-80">{RASTER_REENCODE_ADVICE}</p>
        </div>
    );
}
