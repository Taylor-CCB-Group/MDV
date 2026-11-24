import type { ElementIssue, StoreIssue } from "@spatialdata/core";
import type { EntryStatus } from "@spatialdata/vis";

import { spatialDataUrlForRegion } from "./shared_spatial_data";

/**
 * Data quality: what sd.js reports about how a project's SpatialData stores are written (no
 * pyramid, oversized chunks, no length unit), as opposed to anything about how a layer is set
 * up or loaded.
 *
 * Settings › Data quality lists them for every store. Elsewhere they stay quiet: the layer
 * dialog notes them on one collapsed line, or not at all (`dataQualityNotesSetting`), and they
 * never light a chart's menu pip or row icon. `splitStorageIssues` takes them out of a layer's
 * status for that.
 */

/** One issue of an element's storage; a missing length unit is per coordinate system. */
export function elementIssueId(issue: ElementIssue): string {
    return issue.kind === "no-length-unit" ? `${issue.kind}:${issue.coordinateSystem}` : issue.kind;
}

/**
 * A layer's status without its element's storage issues (`element-issues` notices), and those
 * issues. Undefined when nothing else is left; errors and other notices always stay.
 */
export function splitStorageIssues(status: EntryStatus | undefined): {
    status: EntryStatus | undefined;
    storageIssues: ElementIssue[];
} {
    if (!status) return { status: undefined, storageIssues: [] };
    const storageIssues = status.notices.flatMap((notice) => (notice.kind === "element-issues" ? notice.issues : []));
    if (storageIssues.length === 0) return { status, storageIssues };
    const notices = status.notices.filter((notice) => notice.kind !== "element-issues");
    const left = status.errors.length > 0 || notices.length > 0 ? { errors: status.errors, notices } : undefined;
    return { status: left, storageIssues };
}

/** A store's issues by element path, elements in the order sd.js listed them. */
export function issuesByElement(issues: readonly StoreIssue[]): Map<string, ElementIssue[]> {
    const byPath = new Map<string, ElementIssue[]>();
    for (const { path, issue } of issues) byPath.set(path, [...(byPath.get(path) ?? []), issue]);
    return byPath;
}

type RegionsLike = { all_regions?: Record<string, unknown> } | undefined | null;

/** Every SpatialData store the project's regions point at, once each, with the datasources using it. */
export function projectSpatialStores(
    dataSources: readonly { name: string; dataStore?: { regions?: RegionsLike } }[],
): { url: string; dataSources: string[] }[] {
    const stores = new Map<string, Set<string>>();
    for (const { name, dataStore } of dataSources) {
        for (const region of Object.values(dataStore?.regions?.all_regions ?? {})) {
            const url = spatialDataUrlForRegion(region);
            if (!url) continue;
            stores.set(url, (stores.get(url) ?? new Set()).add(name));
        }
    }
    return [...stores].map(([url, names]) => ({ url, dataSources: [...names] }));
}
