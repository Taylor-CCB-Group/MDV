import { Button, Typography } from "@mui/material";
import type { ElementIssue, SpatialData, StoreIssue } from "@spatialdata/core";
import { observer } from "mobx-react-lite";
import { useEffect, useState, useSyncExternalStore } from "react";

import {
    dataQualityKey,
    elementIssueId,
    isDataQualityDismissed,
    issuesByElement,
    projectSpatialStores,
    setDataQualityDismissed,
} from "@/react/spatialdata/data_quality";
import { describeElementIssue, remediesFor } from "@/react/spatialdata/entry_status";
import { sharedSpatialDataReader } from "@/react/spatialdata/spatial_data_reader";
import { useChartManager } from "../hooks";

/**
 * Settings › Data quality: how each of the project's SpatialData stores is written, as sd.js
 * judges it (`SpatialData.assessIssues`), and what to do about it. Issues dismissed in a layer
 * dialog are listed too, marked, and can be shown there again.
 *
 * Stores are read through the page's shared reader. One a chart has already opened is listed
 * at once; any other waits for "Check", so opening Settings never starts reading a remote store.
 *
 * Each issue carries its `Remedy` as data. Fixing a store from here (asking the server to
 * rewrite it with spatialdata-js-util) would act on that; for now the remedy is described.
 */
export const DataQualitySection = observer(function DataQualitySection() {
    const { dataSources } = useChartManager();
    const stores = projectSpatialStores(dataSources ?? []);
    if (stores.length === 0) {
        return (
            <Typography variant="body2" color="text.secondary">
                This project has no SpatialData stores.
            </Typography>
        );
    }
    return (
        <div className="grid gap-3 pb-2">
            {stores.map((store) => (
                <StoreQuality key={store.url} url={store.url} dataSources={store.dataSources} />
            ))}
        </div>
    );
});

/** Whether the project has any SpatialData store, for showing the section at all. */
export function useHasSpatialStores(): boolean {
    const { dataSources } = useChartManager();
    return projectSpatialStores(dataSources ?? []).length > 0;
}

function storeName(url: string): string {
    const parts = url.replace(/\/+$/, "").split("/");
    return parts.at(-1) || url;
}

function useStoreRead(url: string) {
    const reader = sharedSpatialDataReader();
    return useSyncExternalStore(reader.subscribe, () => reader.peek(url)?.state);
}

function useStoreIssues(spatialData: SpatialData | null) {
    const [result, setResult] = useState<{ spatialData: SpatialData; issues: StoreIssue[] | Error }>();
    useEffect(() => {
        if (!spatialData) return;
        let cancelled = false;
        spatialData.assessIssues().then(
            (issues) => {
                if (!cancelled) setResult({ spatialData, issues: [...issues] });
            },
            (error) => {
                if (!cancelled)
                    setResult({ spatialData, issues: error instanceof Error ? error : new Error(String(error)) });
            },
        );
        return () => {
            cancelled = true;
        };
    }, [spatialData]);
    return result && result.spatialData === spatialData ? result.issues : undefined;
}

function StoreQuality({ url, dataSources }: { url: string; dataSources: string[] }) {
    const state = useStoreRead(url);
    const issues = useStoreIssues(state?.spatialData ?? null);
    const reader = sharedSpatialDataReader();
    return (
        <section className="grid gap-1">
            <div className="flex items-baseline gap-2">
                <Typography variant="subtitle2" title={url}>
                    {storeName(url)}
                </Typography>
                <Typography variant="caption" color="text.secondary">
                    {dataSources.join(", ")}
                </Typography>
            </div>
            {!state ? (
                <div className="flex items-center gap-2">
                    <Typography variant="body2" color="text.secondary">
                        Not opened by any chart yet.
                    </Typography>
                    <Button size="small" variant="outlined" onClick={() => reader.read(url)}>
                        Check
                    </Button>
                </div>
            ) : state.error ? (
                <div className="flex items-center gap-2">
                    <Typography variant="body2" color="error">
                        Could not read the store: {state.error.message}
                    </Typography>
                    <Button size="small" variant="outlined" onClick={() => reader.retry()}>
                        Retry
                    </Button>
                </div>
            ) : state.loading || issues === undefined ? (
                <Typography variant="body2" color="text.secondary">
                    Checking…
                </Typography>
            ) : issues instanceof Error ? (
                <Typography variant="body2" color="error">
                    Could not check the store: {issues.message}
                </Typography>
            ) : issues.length === 0 ? (
                <Typography variant="body2" color="text.secondary">
                    No issues found.
                </Typography>
            ) : (
                <StoreIssueList url={url} issues={issues} />
            )}
        </section>
    );
}

type IssueRow = { id: string; label: string; issue: ElementIssue; paths: string[] };

/**
 * One row per layout issue of an element, and one per coordinate system that has no length
 * unit, naming the elements it affects: that issue is about the coordinate system, and every
 * element with a transform to it reports it.
 */
function issueRows(issues: readonly StoreIssue[]): IssueRow[] {
    const rows: IssueRow[] = [];
    const unitRows = new Map<string, IssueRow>();
    for (const [path, elementIssues] of issuesByElement(issues)) {
        for (const issue of elementIssues) {
            const id = elementIssueId(issue);
            if (issue.kind === "no-length-unit") {
                const row = unitRows.get(id);
                if (row) row.paths.push(path);
                else {
                    const created = {
                        id,
                        label: "Length unit",
                        issue,
                        paths: [path],
                    };
                    unitRows.set(id, created);
                    rows.push(created);
                }
                continue;
            }
            rows.push({ id, label: path, issue, paths: [path] });
        }
    }
    return rows;
}

const StoreIssueList = observer(function StoreIssueList({ url, issues }: { url: string; issues: StoreIssue[] }) {
    const rows = issueRows(issues);
    return (
        <ul className="grid gap-2">
            {rows.map((row) => {
                // dismissed per element: a coordinate system's row can be dismissed for some of its elements
                const dismissedPaths = row.paths.filter((path) =>
                    isDataQualityDismissed(dataQualityKey(url, path, row.id)),
                );
                const all = dismissedPaths.length === row.paths.length;
                const some = dismissedPaths.length > 0;
                const setAll = (value: boolean) => {
                    for (const path of row.paths) setDataQualityDismissed(dataQualityKey(url, path, row.id), value);
                };
                const grouped = row.paths.length > 1 || row.label !== row.paths[0];
                return (
                    <li key={`${row.label}:${row.id}`} className={`grid gap-0.5 ${all ? "opacity-60" : ""}`}>
                        <div className="flex items-start gap-2">
                            <Typography variant="body2" className="grow">
                                <span className="font-medium">{row.label}</span>
                                {all ? " (dismissed)" : some ? ` (dismissed for ${dismissedPaths.join(", ")})` : ""}:{" "}
                                {describeElementIssue(row.issue)}
                            </Typography>
                            {!all && (
                                <Button
                                    size="small"
                                    sx={{ flexShrink: 0 }}
                                    onClick={() => setAll(true)}
                                    title="Stop showing it in the layer dialog and the chart status, in this browser"
                                >
                                    Dismiss
                                </Button>
                            )}
                            {some && (
                                <Button
                                    size="small"
                                    sx={{ flexShrink: 0 }}
                                    onClick={() => setAll(false)}
                                    title="Show it in the layer dialog and the chart status again"
                                >
                                    Show again
                                </Button>
                            )}
                        </div>
                        {grouped && (
                            <Typography variant="caption" color="text.secondary">
                                {row.paths.join(", ")}
                            </Typography>
                        )}
                        <Typography variant="caption" color="text.secondary">
                            {remediesFor([row.issue]).join(" ")}
                        </Typography>
                    </li>
                );
            })}
        </ul>
    );
});
