import { describeElementIssue, type ElementIssue, type EntryNotice, type Remedy } from "@spatialdata/core";
import type { EntryStatus, EntryStatusError } from "@spatialdata/vis";

/**
 * What sd.js reports about each SpatialData layer (`EntryStatus`, through the renderer's
 * `onEntryStatusChange`): why it could not load, and what is degraded but still drawing.
 * MDV keeps these per entry and words the parts sd.js leaves to the host, such as how to fix
 * an element's storage (`Remedy`).
 *
 * Distinct from `layer_failures.ts`, which catches layers that fail while *drawing* and turns
 * them off. A load error here leaves the layer on: sd.js resolves it again when retried.
 */

export type ElementIssuesNotice = Extract<EntryNotice, { kind: "element-issues" }>;

/** How to fix an element's storage, in MDV's words. */
export function describeRemedy(remedy: Remedy): string {
    switch (remedy.op) {
        case "build-pyramid":
            return `Re-encode with a multiscale pyramid: downsample by ${remedy.factor} until the coarsest level is about ${remedy.coarsestAtMostPixels} px or smaller.`;
        case "rechunk":
            return `Rechunk to ${remedy.chunkPixels.min}–${remedy.chunkPixels.max} px chunks (at most ${Math.round(remedy.decodedBytesAtMost / (1024 * 1024))} MB decoded each).`;
        case "declare-length-unit":
            return `Give the spatial axes of coordinate system '${remedy.coordinateSystem}' a length unit, such as micrometer.`;
    }
}

/** One line per distinct remedy, in the order their issues came. */
export function remediesFor(issues: readonly ElementIssue[]): string[] {
    return [...new Set(issues.map((issue) => describeRemedy(issue.remedy)))];
}

export { describeElementIssue };

/**
 * A status as MDV shows it for a layer in `coordinateSystem`.
 *
 * - An element's issues include a missing length unit for every coordinate system it has a
 *   transform to; only the one the chart shows matters here.
 * - Points' notices are left out, as sd.js's own viewer does: the points panel reports their
 *   truncation itself.
 *
 * Undefined when nothing is left to report.
 */
export function entryStatusForDisplay(
    status: EntryStatus | undefined,
    coordinateSystem: string | null,
    elementType?: string,
): EntryStatus | undefined {
    if (!status) return undefined;
    const notices =
        elementType === "points"
            ? []
            : status.notices.flatMap((notice): EntryNotice[] => {
                  if (notice.kind !== "element-issues") return [notice];
                  const issues = notice.issues.filter(
                      (issue) => issue.kind !== "no-length-unit" || issue.coordinateSystem === coordinateSystem,
                  );
                  if (issues.length === 0) return [];
                  if (issues.length === notice.issues.length) return [notice];
                  return [{ ...notice, issues, message: issues.map(describeElementIssue).join(" ") }];
              });
    if (status.errors.length === 0 && notices.length === 0) return undefined;
    return { errors: status.errors, notices };
}

/** Identifies a load error for dismissal: the same failure reported again stays dismissed. */
export function entryStatusErrorKey(entryId: string, error: EntryStatusError): string {
    return `${entryId}\0${error.resource}\0${error.error.kind}\0${error.error.message}`;
}

/** Whether a status says the element is stored in a way that makes it slow. */
export function hasLayoutIssues(status: EntryStatus): boolean {
    return status.notices.some(
        (notice) => notice.kind === "element-issues" && notice.issues.some((issue) => issue.kind !== "no-length-unit"),
    );
}
