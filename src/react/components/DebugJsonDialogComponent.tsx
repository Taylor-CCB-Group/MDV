//import BaseChart from "@/charts/BaseChart";
import JsonView from "react18-json-view";
import "react18-json-view/src/style.css";
// If dark mode is needed, import `dark.css`.
import "react18-json-view/src/dark.css";
import { observer } from "mobx-react-lite";
import { type ReactNode, useEffect, useMemo, useState } from "react";
import { useDebounce } from "use-debounce";
import "../../utilities/css/JsonDialogStyles.css";
import ClearIcon from "@mui/icons-material/Clear";
import { Box, Button, Divider, IconButton, InputAdornment, Link, Paper, TextField, Typography } from "@mui/material";
import type BaseChart from "@/charts/BaseChart";
import { Accordion, AccordionContent, AccordionItem, AccordionTrigger } from "@/components/ui/accordion";
import { browserSettings } from "@/lib/browserSettings";
import type ValidationFindingsStore from "@/lib/ValidationFindingsStore";
import type { ValidationFinding } from "@/lib/ValidationFindingsStore";
import { DEBUG_JSON_REPORT_MESSAGE, MDV_EMAIL } from "@/utilities/constants";
import { ChartProvider } from "../context";
import { useChartManager } from "../hooks";
import BrowserSettingsSection, { browserSettingMatchesQuery } from "./BrowserSettingsSection";
import CustomTooltip from "./CustomTooltip";
import { filterJSON, filterJsonSection, type JsonSection, matchesQuery, splitJsonSections } from "./debugDialogSearch";
import { AbstractComponent, collectDisposers } from "./SettingsDialogComponent";

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
                <Typography sx={{ fontWeight: 700 }}>{finding.subject}</Typography>
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
            {latest?.rawConfig !== undefined && <FindingDetailsJson label="Raw config" value={latest.rawConfig} />}
            {latest?.details !== undefined && <FindingDetailsJson label="Details" value={latest.details} />}
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
                            {subject}: {findings.length} grouped {findings.length === 1 ? "finding" : "findings"},{" "}
                            {repeatCount} total {repeatCount === 1 ? "report" : "reports"}
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

export function ValidationIssuesSection({ validationFindings }: { validationFindings?: ValidationFindingsStore }) {
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
                <Button variant="outlined" size="small" onClick={() => validationFindings.clearAll()}>
                    Clear
                </Button>
            </Box>
            <ValidationFindingGroup title="Charts" findingsBySubject={snapshot.charts} />
            <ValidationFindingGroup title="Datasources" findingsBySubject={snapshot.datasources} />
        </>
    );
}

const PREFERENCES_TITLE = "Preferences";
const DEBUG_GROUP_TITLE = "Debug and diagnostics";
const TECHNICAL_SETTINGS_TITLE = "Technical settings";

/**
 * An accordion section of the dialog, styled like the folders in the chart settings dialog. While a search
 * is active it opens by default (so matches are visible); the user can still toggle it, and those choices
 * are reset when the search changes. Content stays mounted while closed.
 */
function DebugSection({
    id,
    title,
    defaultExpanded = false,
    searching,
    expandedById,
    setExpanded,
    children,
}: {
    id: string;
    title: string;
    defaultExpanded?: boolean;
    children: ReactNode;
} & SectionControlProps) {
    const expanded = expandedById[id] ?? (searching || defaultExpanded);
    return (
        <Paper
            variant="outlined"
            className="w-full"
            sx={{
                backgroundColor: (theme) =>
                    theme.palette.mode === "dark" ? "rgba(255,255,255,0.03)" : "rgba(0,0,0,0.02)",
                borderColor: (theme) => theme.palette.divider,
                borderWidth: 1,
                borderStyle: "solid",
                borderRadius: 1,
                overflow: "hidden",
                paddingX: 2,
                paddingBottom: 0,
                mb: 1,
            }}
        >
            <Accordion
                type="single"
                collapsible
                className="w-full"
                value={expanded ? id : ""}
                onValueChange={(value) => setExpanded(id, value === id)}
            >
                <AccordionItem value={id} className="border-b-0">
                    <AccordionTrigger>{title}</AccordionTrigger>
                    {/* forceMount keeps the content (and its state) mounted while closed; Radix doesn't hide it then, so do it here. */}
                    <AccordionContent forceMount hidden={!expanded}>
                        {children}
                    </AccordionContent>
                </AccordionItem>
            </Accordion>
        </Paper>
    );
}

type SectionControlProps = {
    searching: boolean;
    expandedById: Record<string, boolean>;
    setExpanded: (id: string, expanded: boolean) => void;
};

/** A JSON section with its own filter box; a group renders its children as nested sections. */
function JsonSectionView({
    section,
    idPrefix,
    ...sectionProps
}: { section: JsonSection; idPrefix: string } & SectionControlProps) {
    const id = `${idPrefix}${section.title}`;
    if (section.children) {
        return (
            <DebugSection id={id} title={section.title} {...sectionProps}>
                {section.children.map((child) => (
                    <JsonSectionView key={child.title} section={child} idPrefix={`${id}/`} {...sectionProps} />
                ))}
            </DebugSection>
        );
    }
    return (
        <DebugSection id={id} title={section.title} {...sectionProps}>
            <FilteredJson json={section.json} />
        </DebugSection>
    );
}

function FilteredJson({ json }: { json: unknown }) {
    const [filter, setFilter] = useState("");
    const [debouncedFilter] = useDebounce(filter, 300);
    const filtered = useMemo(
        () => (typeof json === "object" && json !== null ? filterJSON(json, debouncedFilter) : json),
        [json, debouncedFilter],
    );
    return (
        <>
            <TextField
                fullWidth
                size="small"
                variant="outlined"
                label="Filter"
                value={filter}
                onChange={(e) => setFilter(e.target.value)}
                sx={{ mb: 1 }}
            />
            <JsonView
                src={filtered}
                style={{
                    fontSize: "0.875rem",
                    fontFamily: "ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace",
                    padding: "10px",
                }}
                collapsed={1}
            />
        </>
    );
}

const DebugJsonDialogComponent = observer(function DebugJsonDialogComponent({
    json,
    header,
    showValidationSection = false,
    chartPanel,
    groupKeys,
    jsonTitle,
    chart,
}: {
    json: any;
    header?: string;
    showValidationSection?: boolean;
    /**
     * Show `json` whole, as one JSON view under this title, open and just below the chart status, instead of
     * splitting it into a section per top-level key under Debug and diagnostics. For a chart's config.
     */
    jsonTitle?: string;
    /** Top-level keys of `json` to show as groups of sections (e.g. `project: { datasources, views, state }`). */
    groupKeys?: readonly string[];
    /** The chart's own status and problems (`BaseChart.getDebugPanel`), shown first. */
    chartPanel?: ReactNode;
    /**
     * Set for a chart's dialog. That shows the chart's status, its technical settings (`getDebugSettings`) and
     * its config; not the project-wide Preferences or Debug and diagnostics.
     */
    chart?: BaseChart<any>;
}) {
    const [search, setSearch] = useState("");
    const [debouncedSearch] = useDebounce(search, 300); //why is this returning a tuple?
    const query = debouncedSearch.trim();
    const searching = query.length > 0;
    // The user's own open/close choices; reset whenever the search changes, so matches open by default.
    const [expandedState, setExpandedState] = useState<{ query: string; byId: Record<string, boolean> }>({
        query: "",
        byId: {},
    });
    const expandedById = expandedState.query === query ? expandedState.byId : {};
    const setExpanded = (id: string, expanded: boolean) =>
        setExpandedState({ query, byId: { ...expandedById, [id]: expanded } });

    // A search matching the group's own name shows everything in it.
    const debugQuery = matchesQuery(query, DEBUG_GROUP_TITLE) ? "" : query;

    const { validationFindings } = useChartManager();
    const jsonSections = useMemo(
        () => (jsonTitle ? [] : splitJsonSections(omitValidationFindings(json), groupKeys)),
        [json, groupKeys, jsonTitle],
    );
    const wholeJson = useMemo(
        () => (jsonTitle ? filterJsonSection({ title: jsonTitle, json: omitValidationFindings(json) }, query) : null),
        [json, jsonTitle, query],
    );
    const filteredJsonSections = useMemo(
        () =>
            jsonSections.flatMap((section) => {
                const filtered = filterJsonSection(section, debugQuery);
                return filtered ? [filtered] : [];
            }),
        [jsonSections, debugQuery],
    );
    const filteredBrowserSettings = useMemo(
        () =>
            matchesQuery(query, PREFERENCES_TITLE)
                ? browserSettings
                : browserSettings.filter((setting) => browserSettingMatchesQuery(setting, query)),
        [query],
    );
    const technicalSettings = useMemo(() => chart?.getDebugSettings?.() ?? [], [chart]);
    useEffect(() => {
        const disposers = collectDisposers(technicalSettings);
        return () => {
            for (const dispose of disposers) dispose();
        };
    }, [technicalSettings]);
    const filteredTechnicalSettings = useMemo(
        () =>
            matchesQuery(query, TECHNICAL_SETTINGS_TITLE)
                ? technicalSettings
                : technicalSettings.filter((setting) => matchesQuery(query, setting.label)),
        [technicalSettings, query],
    );
    // If there are validation issues when the dialog opens, lead with them instead of Preferences. Captured once,
    // so that clearing the findings doesn't rearrange the dialog under the user.
    const [showFindingsFirst] = useState(() => showValidationSection && Boolean(validationFindings?.hasAny));
    const injectFaultAndReload = (mode: DataLoaderFaultMode) => {
        if (setDataLoaderFaultMode(mode)) {
            window.location.reload();
        }
    };

    const sectionProps = { searching, expandedById, setExpanded };
    const showChartPanel = Boolean(chartPanel) && matchesQuery(query, "Chart status");
    const showReport = matchesQuery(debugQuery, "Report an issue or a bug", MDV_EMAIL, ...DEBUG_JSON_REPORT_MESSAGE);
    const showFaultInjection = matchesQuery(
        debugQuery,
        "Fault Injection",
        "get_data return empty data corrupt data timeout",
    );
    const showValidation =
        showValidationSection && matchesQuery(debugQuery, "Validation issues", "validation findings");
    // a chart's dialog is about that chart; the project-wide sections are in the main Settings dialog
    const showBrowserSettings = !chart && filteredBrowserSettings.length > 0;
    const showTechnicalSettings = filteredTechnicalSettings.length > 0;
    const showDebugGroup =
        !chart && (showReport || showFaultInjection || showValidation || filteredJsonSections.length > 0);
    const nothingMatches =
        searching && !showChartPanel && !showTechnicalSettings && !wholeJson && !showBrowserSettings && !showDebugGroup;

    return (
        <div className="w-full max-h-[80vh] overflow-auto" style={{ userSelect: "text" }}>
            {header && <h2 className="text-lg font-semibold mb-4">{header}</h2>}
            <Box sx={{ p: 1 }}>
                <TextField
                    autoFocus
                    fullWidth
                    size="small"
                    variant="outlined"
                    label="Search Settings, Sections and Data"
                    value={search}
                    onChange={(e) => setSearch(e.target.value)}
                    sx={{ backgroundColor: "transparent" }}
                    slotProps={{
                        input: {
                            endAdornment:
                                search.length > 0 ? (
                                    <InputAdornment position="end">
                                        <IconButton
                                            aria-label="clear search"
                                            onMouseDown={(e) => e.preventDefault()}
                                            onClick={() => setSearch("")}
                                        >
                                            <ClearIcon />
                                        </IconButton>
                                    </InputAdornment>
                                ) : null,
                        },
                    }}
                />
            </Box>
            {nothingMatches && (
                <Typography variant="body1" color="text.secondary" sx={{ marginTop: 1, p: 1 }}>
                    Nothing matches your search.
                </Typography>
            )}
            <Divider sx={{ marginY: 0.75, opacity: 0.6 }} />
            {showChartPanel && (
                <DebugSection id="chart-status" title="Chart status" defaultExpanded {...sectionProps}>
                    {chartPanel}
                </DebugSection>
            )}
            {chart && showTechnicalSettings && (
                <DebugSection
                    id="technical-settings"
                    title={TECHNICAL_SETTINGS_TITLE}
                    defaultExpanded
                    {...sectionProps}
                >
                    <ChartProvider chart={chart}>
                        {filteredTechnicalSettings.map((setting) => (
                            <AbstractComponent key={`${setting.type}:${setting.label}`} props={setting} />
                        ))}
                    </ChartProvider>
                </DebugSection>
            )}
            {wholeJson && (
                <DebugSection id="json" title={wholeJson.title} defaultExpanded {...sectionProps}>
                    <FilteredJson json={wholeJson.json} />
                </DebugSection>
            )}
            {showBrowserSettings && (
                <DebugSection
                    id="preferences"
                    title={PREFERENCES_TITLE}
                    defaultExpanded={!showFindingsFirst}
                    {...sectionProps}
                >
                    <BrowserSettingsSection settings={filteredBrowserSettings} />
                </DebugSection>
            )}
            {showDebugGroup && (
                <DebugSection
                    id="debug"
                    title={DEBUG_GROUP_TITLE}
                    defaultExpanded={showFindingsFirst}
                    {...sectionProps}
                >
                    {showReport && (
                        <DebugSection id="report" title="Report an issue or a bug" {...sectionProps}>
                            <Typography>
                                {DEBUG_JSON_REPORT_MESSAGE[0]}{" "}
                                <Link
                                    href={`mailto:${MDV_EMAIL}?subject=MDV%20Bug%20Report`}
                                    sx={{
                                        textDecoration: "none",
                                        "&.MuiLink-root": { color: "info.main" },
                                    }}
                                >
                                    {MDV_EMAIL}
                                </Link>
                            </Typography>
                            <Typography sx={{ mt: 1 }}>{DEBUG_JSON_REPORT_MESSAGE[1]}</Typography>
                            <Typography sx={{ mt: 1 }}>{DEBUG_JSON_REPORT_MESSAGE[2]}</Typography>
                        </DebugSection>
                    )}
                    {showFaultInjection && (
                        <DebugSection id="fault-injection" title="Fault Injection" {...sectionProps}>
                            <Typography sx={{ mb: 1 }}>
                                Select a one-shot fault to apply to the next <code>/get_data</code> request. The project
                                reloads automatically.
                            </Typography>
                            <Box sx={{ display: "flex", flexWrap: "wrap", gap: 1 }}>
                                <CustomTooltip tooltipText="The frontend discards one loaded column batch and returns no data.">
                                    <Button
                                        variant="outlined"
                                        size="small"
                                        onClick={() => injectFaultAndReload("empty")}
                                    >
                                        Return empty data
                                    </Button>
                                </CustomTooltip>
                                <CustomTooltip tooltipText="The frontend replaces one response body with invalid bytes.">
                                    <Button
                                        variant="outlined"
                                        size="small"
                                        onClick={() => injectFaultAndReload("corrupt")}
                                    >
                                        Corrupt data
                                    </Button>
                                </CustomTooltip>
                                <CustomTooltip tooltipText="Client-side one-shot delay, then failure.">
                                    <Button
                                        variant="outlined"
                                        size="small"
                                        onClick={() => injectFaultAndReload("timeout")}
                                    >
                                        Timeout
                                    </Button>
                                </CustomTooltip>
                            </Box>
                        </DebugSection>
                    )}
                    {showValidation && (
                        <DebugSection
                            id="validation-issues"
                            title="Validation issues"
                            defaultExpanded={showFindingsFirst}
                            {...sectionProps}
                        >
                            <ValidationIssuesSection validationFindings={validationFindings} />
                        </DebugSection>
                    )}
                    {filteredJsonSections.map((section) => (
                        <JsonSectionView key={section.title} section={section} idPrefix="json:" {...sectionProps} />
                    ))}
                </DebugSection>
            )}
        </div>
    );
});

export default DebugJsonDialogComponent;
