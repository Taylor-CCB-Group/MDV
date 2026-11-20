import type { ElementIssue, EntryNotice, StoreIssue } from "@spatialdata/core";
import type { RenderStackSpatialElementType } from "@spatialdata/layers";
import type { EntryStatus } from "@spatialdata/vis";
import { action, observable } from "mobx";

import { spatialDataUrlForRegion } from "./shared_spatial_data";

/**
 * Data quality: what sd.js reports about how a project's SpatialData stores are written, and
 * which of those reports the user has dismissed.
 *
 * A dismissal belongs to this browser (localStorage), like the Preferences settings: it stops
 * the layer dialog and the chart's menu pip from repeating something the user knows about,
 * and is never saved with a project. Settings › Data quality still lists every issue, marking
 * the dismissed ones, and can restore them. A dismissal is keyed by the store, the element
 * and the issue (`dataQualityKey`), so it holds for every chart and layer showing that
 * element, and a different issue on the same element still shows.
 */

const STORAGE_KEY = "mdv.dataQuality.dismissed";

const ELEMENT_FOLDERS: Record<RenderStackSpatialElementType, string> = {
    image: "images",
    labels: "labels",
    points: "points",
    shapes: "shapes",
};

/** An element's path in its store, as sd.js's `StoreIssue.path` names it (`images/he`). */
export function elementPath(elementType: RenderStackSpatialElementType, elementKey: string): string {
    return `${ELEMENT_FOLDERS[elementType]}/${elementKey}`;
}

/** One issue of an element's storage; a missing length unit is per coordinate system. */
export function elementIssueId(issue: ElementIssue): string {
    return issue.kind === "no-length-unit" ? `${issue.kind}:${issue.coordinateSystem}` : issue.kind;
}

/** What a dismissal is keyed by: the store, the element, and the issue or notice kind. */
export function dataQualityKey(storeUrl: string, path: string, id: string): string {
    return JSON.stringify([storeUrl, path, id]);
}

function readDismissed(): string[] {
    try {
        const parsed: unknown = JSON.parse(window.localStorage?.getItem(STORAGE_KEY) ?? "[]");
        return Array.isArray(parsed) ? parsed.filter((item): item is string => typeof item === "string") : [];
    } catch {
        return [];
    }
}

const dismissed = observable.set<string>(typeof window === "undefined" ? [] : readDismissed());

export function isDataQualityDismissed(key: string): boolean {
    return dismissed.has(key);
}

export const setDataQualityDismissed = action((key: string, value: boolean) => {
    if (value) dismissed.add(key);
    else dismissed.delete(key);
    try {
        window.localStorage?.setItem(STORAGE_KEY, JSON.stringify([...dismissed]));
    } catch (error) {
        console.warn("Could not save dismissed data quality issues", error);
    }
});

/** Each dismissable part of a notice: an element's issues one by one, any other notice whole. */
export function noticeParts(notice: EntryNotice): string[] {
    return notice.kind === "element-issues" ? notice.issues.map(elementIssueId) : [notice.kind];
}

/**
 * A layer's status without what has been dismissed for its element: notices, and the
 * dismissed issues of an `element-issues` notice. Errors always stay. Undefined when nothing
 * is left; `dismissed` counts what was taken out.
 */
export function withoutDismissed(
    status: EntryStatus | undefined,
    storeUrl: string | null,
    path: string | null,
): { status: EntryStatus | undefined; dismissed: number } {
    if (!status) return { status: undefined, dismissed: 0 };
    if (!storeUrl || !path) return { status, dismissed: 0 };
    const isDismissed = (id: string) => dismissed.has(dataQualityKey(storeUrl, path, id));
    let count = 0;
    const notices = status.notices.flatMap((notice): EntryNotice[] => {
        if (notice.kind !== "element-issues") {
            if (!isDismissed(notice.kind)) return [notice];
            count++;
            return [];
        }
        const issues = notice.issues.filter((issue) => !isDismissed(elementIssueId(issue)));
        count += notice.issues.length - issues.length;
        if (issues.length === 0) return [];
        return [issues.length === notice.issues.length ? notice : { ...notice, issues }];
    });
    const left = status.errors.length > 0 || notices.length > 0 ? { errors: status.errors, notices } : undefined;
    return { status: left, dismissed: count };
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
