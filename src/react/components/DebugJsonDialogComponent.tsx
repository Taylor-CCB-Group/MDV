//import BaseChart from "@/charts/BaseChart";
import JsonView from "react18-json-view";
import "react18-json-view/src/style.css";
// If dark mode is needed, import `dark.css`.
import "react18-json-view/src/dark.css";
import { observer } from "mobx-react-lite";
import { useMemo, useState } from "react";
import { useDebounce } from "use-debounce";
import { vivLoaderCacheTelemetryObservable } from "../viv_loader_cache";
import "../../utilities/css/JsonDialogStyles.css";
import { Accordion, AccordionDetails, AccordionSummary, Box, Button, Divider, Link, Typography } from "@mui/material";
import ExpandMoreIcon from '@mui/icons-material/ExpandMore';
import { DEBUG_JSON_REPORT_MESSAGE, MDV_EMAIL } from "@/utilities/constants";
import { useChartManager } from "../hooks";
import CustomTooltip from "./CustomTooltip";
import type ValidationFindingsStore from "@/lib/ValidationFindingsStore";
import type { ValidationFinding } from "@/lib/ValidationFindingsStore";

type JSONObject = { [key: string]: any };
type DataLoaderFaultMode = "empty" | "corrupt" | "timeout";

const DATA_LOADER_FAULT_STORAGE_KEY = "mdvDataLoaderFault";

function setDataLoaderFaultMode(mode: DataLoaderFaultMode | null) {
    try {
        if (!mode) {
            window.sessionStorage.removeItem(DATA_LOADER_FAULT_STORAGE_KEY);
            return true;
        }
        window.sessionStorage.setItem(
            DATA_LOADER_FAULT_STORAGE_KEY,
            JSON.stringify({
                mode,
                delayMs: 10_000,
                createdAt: Date.now(),
            }),
        );
        return true;
    } catch (error) {
        console.warn("Failed to write data loader fault config", error);
        return false;
    }
}

/**
 * given a filter string and a json object, return a new json object that is a subset of the input json
 * matching the filter should mean that a given node is included in the output
 * if all (case insensitive) elements of filterArray appear somewhere in the node or its path.
 * credit ChatGPT for the code (copilot failed) 
 * Adendum: the ChatGPT code did not meet the spec, so I had to fix it.
 * This had caused me to spend a lot longer debugging when I was trying to use it before I realised this was at fault.
 * */
function filterJSON(obj: JSONObject, filter: string): JSONObject {
    if (!filter) return obj;
    const filterArray = filter.toLowerCase().split(" ");

    function matchesFilter(value: string, path: string[]): boolean {
        const combinedString = path.concat(value).join(" ").toLowerCase();
        return filterArray.every((f) => combinedString.includes(f));
    }

    function recursiveFilter(
        current: JSONObject,
        path: string[],
    ): JSONObject | null {
        // we should return the whole subtree if the current node matches the filter
        if (matchesFilter("", path)) {
            return current;
        }
        if (typeof current === "string") {
            return matchesFilter(current, path) ? current : null;
        }
        //`typeof current === "object"` probably no longer needed
        if (typeof current === "object" && current !== null) {
            const filteredObj: JSONObject = {};
            let childMatch = false;

            for (const key in current) {
                const result = recursiveFilter(current[key], path.concat(key));
                if (result !== null) {
                    filteredObj[key] = result;
                    childMatch = true;
                }
            }
            if (childMatch || matchesFilter("", path)) {
                return filteredObj;
            }
        }
        return null;
    }

    return recursiveFilter(obj, []) || {};
}

function omitValidationFindings(json: any) {
    if (typeof json !== "object" || json === null || !("validationFindings" in json)) {
        return json;
    }
    const { validationFindings: _validationFindings, ...jsonWithoutValidationFindings } = json;
    return jsonWithoutValidationFindings;
}

function formatSeenAt(timestamp: number) {
    return new Date(timestamp).toISOString();
}

function issuePathLabel(path: string) {
    return path || "(root)";
}

function FindingDetailsJson({ label, value }: { label: string; value: unknown }) {
    return (
        <Box
            component="details"
            sx={{
                mt: 1,
                "& summary": {
                    cursor: "pointer",
                    fontWeight: 600,
                },
            }}
        >
            <summary>{label}</summary>
            <JsonView
                src={value}
                style={{
                    fontSize: "0.8125rem",
                    fontFamily: "ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace",
                    padding: "8px",
                }}
                collapsed={1}
            />
        </Box>
    );
}

function ValidationFindingCard({ finding }: { finding: ValidationFinding }) {
    const latest = finding.latest;

    return (
        <Box
            sx={{
                border: "1px solid var(--border_menu_bar_color)",
                borderRadius: "8px",
                p: 1.5,
                mb: 1,
            }}
        >
            <Box sx={{ display: "flex", flexWrap: "wrap", alignItems: "baseline", gap: 1 }}>
                <Typography sx={{ fontWeight: 700 }}>
                    {finding.subject}
                </Typography>
                <Typography variant="body2" sx={{ opacity: 0.75 }}>
                    Repeated {finding.count} {finding.count === 1 ? "time" : "times"}
                </Typography>
                <Typography variant="body2" sx={{ opacity: 0.75 }}>
                    Latest: {formatSeenAt(finding.lastSeenAt)}
                </Typography>
            </Box>
            <Box sx={{ mt: 1 }}>
                {latest?.issues?.length ? (
                    latest.issues.map((issue, index) => (
                        <Box
                            key={`${issue.path}-${issue.message}-${index}`}
                            sx={{
                                borderLeft: "3px solid #f2c94c",
                                pl: 1,
                                mb: 1,
                            }}
                        >
                            <Typography
                                variant="body2"
                                sx={{
                                    fontFamily: "ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace",
                                    fontWeight: 600,
                                }}
                            >
                                {issuePathLabel(issue.path)}
                                {issue.code ? ` (${issue.code})` : ""}
                            </Typography>
                            <Typography variant="body2">{issue.message}</Typography>
                        </Box>
                    ))
                ) : (
                    <Typography variant="body2">No issue details recorded.</Typography>
                )}
            </Box>
            {latest?.rawConfig !== undefined && (
                <FindingDetailsJson label="Raw config" value={latest.rawConfig} />
            )}
            {latest?.details !== undefined && (
                <FindingDetailsJson label="Details" value={latest.details} />
            )}
        </Box>
    );
}

function ValidationFindingGroup({
    title,
    findingsBySubject,
}: {
    title: string;
    findingsBySubject: Record<string, ValidationFinding[]>;
}) {
    const entries = Object.entries(findingsBySubject).sort(([a], [b]) => a.localeCompare(b));
    if (entries.length === 0) return null;

    const totalFindings = entries.reduce((sum, [, findings]) => sum + findings.length, 0);

    return (
        <Box sx={{ mb: 2 }}>
            <Typography variant="subtitle1" sx={{ fontWeight: 700, mb: 1 }}>
                {title} ({totalFindings})
            </Typography>
            {entries.map(([subject, findings]) => {
                const repeatCount = findings.reduce((sum, finding) => sum + finding.count, 0);
                return (
                    <Box key={subject} sx={{ mb: 1.5 }}>
                        <Typography variant="body2" sx={{ mb: 0.75, opacity: 0.8 }}>
                            {subject}: {findings.length} grouped {findings.length === 1 ? "finding" : "findings"}, {repeatCount} total {repeatCount === 1 ? "report" : "reports"}
                        </Typography>
                        {findings.map((finding) => (
                            <ValidationFindingCard key={finding.signature} finding={finding} />
                        ))}
                    </Box>
                );
            })}
        </Box>
    );
}

export function ValidationIssuesSection({
    validationFindings,
}: {
    validationFindings?: ValidationFindingsStore;
}) {
    const snapshot = validationFindings?.snapshot;
    const hasAnyFindings = Boolean(validationFindings?.hasAny);

    if (!hasAnyFindings || !snapshot || !validationFindings) {
        return <Typography>No validation issues recorded.</Typography>;
    }

    return (
        <>
            <Box sx={{ display: "flex", alignItems: "center", justifyContent: "space-between", gap: 2, mb: 1 }}>
                <Typography>
                    Total: <strong>{snapshot.totalCount}</strong>
                </Typography>
                <Button
                    variant="outlined"
                    size="small"
                    onClick={() => validationFindings.clearAll()}
                >
                    Clear
                </Button>
            </Box>
            <ValidationFindingGroup title="Charts" findingsBySubject={snapshot.charts} />
            <ValidationFindingGroup title="Datasources" findingsBySubject={snapshot.datasources} />
        </>
    );
}

const DebugJsonDialogComponent = observer(function DebugJsonDialogComponent({
    json,
    header,
    showValidationSection = false,
}: {
    json: any;
    header?: string;
    showValidationSection?: boolean;
}) {
    const [filter, setFilter] = useState("");
    const [debouncedFilter] = useDebounce(filter, 300); //why is this returning a tuple?
    const { validationFindings } = useChartManager();
    const jsonWithoutValidationFindings = useMemo(
        () => omitValidationFindings(json),
        [json],
    );
    const jsonWithVivTelemetry = useMemo(
        () => ({
            ...jsonWithoutValidationFindings,
            vivLoaderCacheTelemetry: vivLoaderCacheTelemetryObservable.snapshot,
        }),
        [jsonWithoutValidationFindings, vivLoaderCacheTelemetryObservable.snapshot],
    );
    const filteredJson = useMemo(
        () => filterJSON(jsonWithVivTelemetry, debouncedFilter),
        [jsonWithVivTelemetry, debouncedFilter],
    );
    const hasAnyFindings = Boolean(validationFindings?.hasAny);
    const injectFaultAndReload = (mode: DataLoaderFaultMode) => {
        if (setDataLoaderFaultMode(mode)) {
            window.location.reload();
        }
    };

    return (
        <div className="overflow-auto p-4" style={{ userSelect: "text" }}>
            {header && <h2 className="text-lg font-semibold mb-4">{header}</h2>}
            <Accordion sx={{
                border: "1px solid var(--border_menu_bar_color)",
                mb: 2,
            }}>
                <AccordionSummary expandIcon={<ExpandMoreIcon />}>
                    <Typography variant="h6">
                        Report an issue or a bug
                    </Typography>
                </AccordionSummary>
                <Divider />
                <AccordionDetails sx={{ mt: 1 }}>
                    <Typography>
                        {DEBUG_JSON_REPORT_MESSAGE[0]}
                        {' '}
                        <Link 
                            href={`mailto:${MDV_EMAIL}?subject=MDV%20Bug%20Report`} 
                            sx={{ 
                                textDecoration: "none",
                                "&.MuiLink-root": { color: "info.main" }
                            }}
                        >
                            {MDV_EMAIL}
                        </Link>
                    </Typography>
                    <Typography sx={{mt: 1}}>
                        {DEBUG_JSON_REPORT_MESSAGE[1]}
                    </Typography>
                    <Typography sx={{mt: 1}}>
                        {DEBUG_JSON_REPORT_MESSAGE[2]}
                    </Typography>
                </AccordionDetails>
            </Accordion>

            <Accordion sx={{
                border: "1px solid var(--border_menu_bar_color)",
                mb: 2,
            }}>
                <AccordionSummary expandIcon={<ExpandMoreIcon />}>
                    <Typography variant="h6">
                        Fault Injection
                    </Typography>
                </AccordionSummary>
                <Divider />
                <AccordionDetails sx={{ mt: 1 }}>
                    <Typography sx={{ mb: 1 }}>
                        Select a one-shot fault to apply to the next <code>/get_data</code> request. The project reloads automatically.
                    </Typography>
                    <Box sx={{ display: "flex", flexWrap: "wrap", gap: 1 }}>
                        <CustomTooltip tooltipText="The frontend discards one loaded column batch and returns no data.">
                            <Button variant="outlined" size="small" onClick={() => injectFaultAndReload("empty")}>
                                Return empty data
                            </Button>
                        </CustomTooltip>
                        <CustomTooltip tooltipText="The frontend replaces one response body with invalid bytes.">
                            <Button variant="outlined" size="small" onClick={() => injectFaultAndReload("corrupt")}>
                                Corrupt data
                            </Button>
                        </CustomTooltip>
                        <CustomTooltip tooltipText="Client-side one-shot delay, then failure.">
                            <Button variant="outlined" size="small" onClick={() => injectFaultAndReload("timeout")}>
                                Timeout
                            </Button>
                        </CustomTooltip>
                    </Box>
                </AccordionDetails>
            </Accordion>

            {showValidationSection && (
                <Accordion
                    defaultExpanded={hasAnyFindings}
                    sx={{
                        border: "1px solid var(--border_menu_bar_color)",
                        mb: 2,
                    }}
                >
                    <AccordionSummary expandIcon={<ExpandMoreIcon />}>
                        <Typography variant="h6">
                            Validation issues
                        </Typography>
                    </AccordionSummary>
                    <Divider />
                    <AccordionDetails sx={{ mt: 1 }}>
                        <ValidationIssuesSection validationFindings={validationFindings} />
                    </AccordionDetails>
                </Accordion>
            )}

            <input
                type="text"
                value={filter}
                onChange={(e) => setFilter(e.target.value)}
                placeholder="Filter..."
                className="w-full mb-4 mt-4 p-2 border rounded-md dark:bg-gray-700 dark:border-gray-600"
            />
            <JsonView 
                src={filteredJson} 
                style={{
                    fontSize: '0.875rem',
                    fontFamily: 'ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace',
                    padding: "10px"
                }}
                collapsed={1}
            />
        </div>
    );
});

export default DebugJsonDialogComponent;
