import {
    Alert,
    Button,
    Chip,
    Stack,
    Table,
    TableBody,
    TableCell,
    TableHead,
    TableRow,
    TextField,
    Typography,
} from "@mui/material";
import { useMutation, useQuery, useQueryClient } from "@tanstack/react-query";
import { useId, useState } from "react";
import {
    type JobFormValues,
    type JobsDataSource,
    initialValues,
    paramOptions,
    setParam,
    toSubmitParams,
} from "./jobForm";
import {
    type JobParamSpec,
    type JobParamValue,
    type JobRecord,
    type JobStatus,
    type JobTool,
    fetchExecutors,
    fetchTools,
    submitJob,
} from "./jobsApi";
import { jobsQueryOptions } from "./useJobColumnSync";

const FAILED: JobStatus[] = ["failed", "lost", "stale"];

const EXECUTOR_LABELS: Record<string, string> = { local: "Local", slurm: "Slurm" };

function executorLabel(name: string) {
    return EXECUTOR_LABELS[name] ?? name;
}

type ParamFieldProps = {
    param: JobParamSpec;
    values: JobFormValues;
    sources: JobsDataSource[];
    onChange: (name: string, value: string) => void;
};

function ParamField({ param, values, sources, onChange }: ParamFieldProps) {
    const id = useId();
    const options = paramOptions(param, values, sources);
    const value = values[param.name] ?? "";
    const numeric = param.type === "int" || param.type === "float";

    if (options === null) {
        return (
            <TextField
                id={id}
                label={param.label}
                value={value}
                onChange={(e) => onChange(param.name, e.target.value)}
                size="small"
                fullWidth
                slotProps={{ htmlInput: numeric ? { inputMode: "decimal" } : undefined }}
            />
        );
    }
    return (
        <TextField
            id={id}
            select
            label={param.label}
            value={value}
            onChange={(e) => onChange(param.name, e.target.value)}
            size="small"
            fullWidth
            helperText={options.length === 0 ? "Nothing to choose from on this datasource" : undefined}
            slotProps={{ select: { native: true }, inputLabel: { shrink: true } }}
        >
            {param.type === "column" && (
                <option value="" disabled>
                    Choose a column
                </option>
            )}
            {options.map((o) => (
                <option key={o.value} value={o.value}>
                    {o.label}
                </option>
            ))}
        </TextField>
    );
}

type SelectionFieldProps = {
    param: JobParamSpec;
    source: JobsDataSource;
    selecting: boolean;
    name: string;
    onChoose: (selecting: boolean) => void;
    onName: (name: string) => void;
};

/** Run on the datasource's current selection, which needs a name, or on all of its rows (ADR-0017). */
function SelectionField({ param, source, selecting, name, onChoose, onName }: SelectionFieldProps) {
    const choiceId = useId();
    const nameId = useId();
    const size = source.size.toLocaleString();
    return (
        <>
            <TextField
                id={choiceId}
                select
                label={param.label}
                value={selecting ? "selection" : "all"}
                onChange={(e) => onChoose(e.target.value === "selection")}
                size="small"
                fullWidth
                slotProps={{ select: { native: true }, inputLabel: { shrink: true } }}
            >
                <option value="selection">{`Current selection (${source.filterSize.toLocaleString()} of ${size})`}</option>
                <option value="all">{`All (${size})`}</option>
            </TextField>
            {selecting && (
                <TextField
                    id={nameId}
                    label="Selection name"
                    value={name}
                    onChange={(e) => onName(e.target.value)}
                    size="small"
                    fullWidth
                    helperText={name.trim() ? undefined : "Name the selection to run on it"}
                />
            )}
        </>
    );
}

/** "T cells (2,314)" for a job that ran on a selection, "All" otherwise. */
function cellsLabel(job: JobRecord, tools: JobTool[]) {
    if (job.selection_rows === null) return "All";
    const param = tools.find((t) => t.id === job.tool_id)?.params.find((p) => p.type === "selection");
    const name = param ? job.params[param.name] : undefined;
    return `${name ?? "Selection"} (${job.selection_rows.toLocaleString()})`;
}

function statusColor(status: JobStatus) {
    if (status === "done") return "success";
    if (FAILED.includes(status)) return "error";
    return "default";
}

type JobsPanelProps = {
    /** project root from `useProject()`: "" in single-project mode, "/project/<id>" otherwise */
    root: string;
    sources: JobsDataSource[];
};

/** Pick a tool, fill in its params, run it, and follow the project's jobs. */
export default function JobsPanel({ root, sources }: JobsPanelProps) {
    const queryClient = useQueryClient();
    const toolsQuery = useQuery({ queryKey: ["jobs-tools", root], queryFn: () => fetchTools(root) });
    const jobsQuery = useQuery(jobsQueryOptions(root));
    const executorsQuery = useQuery({ queryKey: ["jobs-executors", root], queryFn: () => fetchExecutors(root) });
    const submit = useMutation({
        mutationFn: async ({
            toolId,
            params,
            executor,
            selectedFrom,
        }: {
            toolId: string;
            params: Record<string, JobParamValue>;
            executor?: string;
            /** the datasource whose current selection the job runs on */
            selectedFrom?: JobsDataSource;
        }) => {
            const rows = selectedFrom ? Array.from(await selectedFrom.getFilteredIndices()) : undefined;
            return submitJob(root, toolId, params, executor, rows);
        },
        onSuccess: () => queryClient.invalidateQueries({ queryKey: ["jobs", root] }),
    });

    const tools = toolsQuery.data ?? [];
    const toolIdSelectId = useId();
    // form values belong to one tool; choosing another tool starts again from its initial values
    const [form, setForm] = useState<{ toolId: string; values: JobFormValues } | null>(null);
    const tool = tools.find((t) => t.id === form?.toolId) ?? tools[0];
    const values = tool && form?.toolId === tool.id ? form.values : tool ? initialValues(tool, sources) : {};
    const [formError, setFormError] = useState<string | null>(null);

    // with a filter active on the job's datasource, the job runs on the selection unless the user
    // chooses all rows (ADR-0017); the dialog is modal, so the filter cannot change while it is open
    const selectionParam = tool?.params.find((p) => p.type === "selection");
    const selectionSource = selectionParam
        ? sources.find((s) => s.name === values[selectionParam.options_from ?? ""])
        : undefined;
    const filtered = selectionSource !== undefined && selectionSource.filterSize < selectionSource.size;
    const [runOnAll, setRunOnAll] = useState(false);
    const selecting = filtered && !runOnAll;
    const needsName = selecting && selectionParam !== undefined && !(values[selectionParam.name] ?? "").trim();

    // the first executor is the server's default (ADR-0016); with only one there is nothing to choose
    const executors = executorsQuery.data ?? [];
    const chooseExecutor = executors.length > 1;
    const executorSelectId = useId();
    const [chosenExecutor, setChosenExecutor] = useState<string | null>(null);
    const executor = chosenExecutor ?? executors[0];

    const toolName = (id: string) => tools.find((t) => t.id === id)?.name ?? id;

    const onRun = () => {
        if (!tool) return;
        setFormError(null);
        submit.reset();
        try {
            submit.mutate({
                toolId: tool.id,
                params: toSubmitParams(tool, values),
                executor: chooseExecutor ? executor : undefined,
                selectedFrom: selecting ? selectionSource : undefined,
            });
        } catch (e) {
            setFormError(e instanceof Error ? e.message : String(e));
        }
    };

    const error =
        formError ??
        submit.error?.message ??
        toolsQuery.error?.message ??
        executorsQuery.error?.message ??
        jobsQuery.error?.message;
    const jobs = jobsQuery.data ?? [];

    return (
        <Stack spacing={2}>
            {tool && (
                <>
                    <TextField
                        id={toolIdSelectId}
                        select
                        label="Tool"
                        value={tool.id}
                        onChange={(e) => {
                            const next = tools.find((t) => t.id === e.target.value);
                            if (next) setForm({ toolId: next.id, values: initialValues(next, sources) });
                        }}
                        size="small"
                        fullWidth
                        slotProps={{ select: { native: true }, inputLabel: { shrink: true } }}
                    >
                        {tools.map((t) => (
                            <option key={t.id} value={t.id}>
                                {t.name}
                            </option>
                        ))}
                    </TextField>
                    <Typography variant="body2" color="text.secondary">
                        {tool.description}
                    </Typography>
                    {tool.params.map((p) =>
                        p.type === "selection" ? (
                            filtered && (
                                <SelectionField
                                    key={`${tool.id}:${p.name}`}
                                    param={p}
                                    source={selectionSource}
                                    selecting={selecting}
                                    name={values[p.name] ?? ""}
                                    onChoose={(choice) => {
                                        setRunOnAll(!choice);
                                        // running on all rows drops the name, and the output name it set
                                        if (!choice) {
                                            setForm({
                                                toolId: tool.id,
                                                values: setParam(tool, values, p.name, "", sources),
                                            });
                                        }
                                    }}
                                    onName={(name) =>
                                        setForm({
                                            toolId: tool.id,
                                            values: setParam(tool, values, p.name, name, sources),
                                        })
                                    }
                                />
                            )
                        ) : (
                            <ParamField
                                key={`${tool.id}:${p.name}`}
                                param={p}
                                values={values}
                                sources={sources}
                                onChange={(name, value) =>
                                    setForm({ toolId: tool.id, values: setParam(tool, values, name, value, sources) })
                                }
                            />
                        ),
                    )}
                    {chooseExecutor && (
                        <TextField
                            id={executorSelectId}
                            select
                            label="Run on"
                            value={executor}
                            onChange={(e) => setChosenExecutor(e.target.value)}
                            size="small"
                            fullWidth
                            slotProps={{ select: { native: true }, inputLabel: { shrink: true } }}
                        >
                            {executors.map((name) => (
                                <option key={name} value={name}>
                                    {executorLabel(name)}
                                </option>
                            ))}
                        </TextField>
                    )}
                    <Button variant="contained" onClick={onRun} disabled={submit.isPending || needsName}>
                        Run
                    </Button>
                </>
            )}
            {error && <Alert severity="error">{error}</Alert>}
            <Typography variant="subtitle1">Jobs</Typography>
            {jobs.length === 0 ? (
                <Typography variant="body2" color="text.secondary">
                    No jobs yet.
                </Typography>
            ) : (
                <Table size="small">
                    <TableHead>
                        <TableRow>
                            <TableCell>Tool</TableCell>
                            <TableCell>Status</TableCell>
                            <TableCell>Cells</TableCell>
                            {chooseExecutor && <TableCell>Run on</TableCell>}
                            <TableCell>Started</TableCell>
                            <TableCell>Details</TableCell>
                        </TableRow>
                    </TableHead>
                    <TableBody>
                        {jobs.map((j) => (
                            <TableRow key={j.job_id}>
                                <TableCell>{toolName(j.tool_id)}</TableCell>
                                <TableCell>
                                    <Chip size="small" label={j.status} color={statusColor(j.status)} />
                                </TableCell>
                                <TableCell>{cellsLabel(j, tools)}</TableCell>
                                {chooseExecutor && (
                                    // a record from before executors had names ran on the default
                                    <TableCell>{executorLabel(j.executor ?? executors[0])}</TableCell>
                                )}
                                <TableCell>{new Date(j.created * 1000).toLocaleString()}</TableCell>
                                <TableCell>{j.error ?? ""}</TableCell>
                            </TableRow>
                        ))}
                    </TableBody>
                </Table>
            )}
        </Stack>
    );
}
