import type { SpatialData } from "@spatialdata/core";
import type { RenderStack } from "@spatialdata/layers";
import { TriangleAlertIcon, XIcon } from "lucide-react";
import { useEffect, useMemo, useState } from "react";
import {
    assessRasterLevels,
    describeRasterIssue,
    loadRasterLevels,
    RASTER_REENCODE_ADVICE,
    type RasterIssue,
} from "../spatialdata/raster_health";

type RasterElementType = "image" | "labels";
type ElementIssues = { elementType: RasterElementType; elementKey: string; issues: RasterIssue[] };
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

function useRasterIssues(spatialData: SpatialData | null | undefined, entriesKey: string) {
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

/**
 * Warns, in the chart itself, when an image or labels element in the render stack is laid
 * out in a way that makes it slow to load and render. See `raster_health.ts`.
 */
export default function RasterDataWarning({
    spatialData,
    entriesKey,
}: {
    spatialData: SpatialData | null | undefined;
    /** From {@link rasterEntriesKey}. */
    entriesKey: string;
}) {
    const issues = useRasterIssues(spatialData, entriesKey);
    const [expanded, setExpanded] = useState(false);
    const [dismissedKey, setDismissedKey] = useState<string | null>(null);
    const issuesKey = useMemo(() => issues.map((item) => `${item.elementType}:${item.elementKey}`).join("|"), [issues]);
    // Dismissal sticks until a different set of elements is flagged.
    if (issues.length === 0 || dismissedKey === issuesKey) return null;

    return (
        <div className="pointer-events-auto absolute bottom-8 right-2 z-[4] flex max-w-[min(28rem,calc(100%-1rem))] flex-col-reverse items-end gap-1 text-xs">
            <div className="flex items-center gap-1 rounded border border-amber-500/60 bg-amber-50/95 text-amber-900 shadow-sm dark:bg-amber-950/90 dark:text-amber-100">
                <button
                    type="button"
                    className="flex flex-1 items-center gap-1.5 px-2 py-1 text-left"
                    aria-expanded={expanded}
                    onClick={() => setExpanded((value) => !value)}
                >
                    <TriangleAlertIcon className="h-3.5 w-3.5 shrink-0" />
                    <span className="font-medium">Slow image data{issues.length > 1 ? ` (${issues.length})` : ""}</span>
                </button>
                <button
                    type="button"
                    className="px-1.5 py-1 opacity-70 hover:opacity-100"
                    aria-label="Dismiss image data warning"
                    onClick={() => setDismissedKey(issuesKey)}
                >
                    <XIcon className="h-3.5 w-3.5" />
                </button>
            </div>
            {expanded && (
                <div className="max-h-[60vh] space-y-2 overflow-y-auto rounded border border-amber-500/60 bg-amber-50/95 px-2 py-1.5 text-amber-900 shadow-sm dark:bg-amber-950/90 dark:text-amber-100">
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
            )}
        </div>
    );
}
