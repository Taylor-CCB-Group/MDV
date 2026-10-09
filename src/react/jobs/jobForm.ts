/**
 * Turns a tool's params (`GET /jobs/tools`) into form state and back into submit params.
 *
 * Form values are kept as strings, the way inputs hold them; `toSubmitParams` converts numeric
 * params when the form is submitted.
 */
import type { JobParamSpec, JobParamType, JobParamValue, JobTool } from "./jobsApi";

/** Param types the user picks from a list rather than typing. */
const PICKERS: JobParamType[] = ["datasource", "subgroup", "column"];

/** The parts of a datasource the jobs form needs, so the form logic does not depend on ChartManager. */
export type JobsDataSource = {
    name: string;
    columns: { field: string; name: string }[];
    /** keys of the `rows_as_columns` subgroups linked from this datasource */
    subgroups: string[];
};

/** The parts of a ChartManager datasource that `toJobsDataSource` reads. */
type AppDataSource = {
    name: string;
    dataStore: {
        columns: { field: string; name: string; subgroup?: string }[];
        links?: Record<string, { rows_as_columns?: { subgroups: Record<string, unknown> } }>;
    };
};

export function toJobsDataSource(ds: AppDataSource): JobsDataSource {
    return {
        name: ds.name,
        // subgroup columns (`gs|CD4|12`) are virtual, so they are not on the server's datasource
        columns: ds.dataStore.columns.filter((c) => !c.subgroup).map((c) => ({ field: c.field, name: c.name })),
        subgroups: Object.values(ds.dataStore.links ?? {}).flatMap((l) =>
            Object.keys(l.rows_as_columns?.subgroups ?? {}),
        ),
    };
}

export type JobFormValues = Record<string, string>;

export type JobParamOption = { value: string; label: string };

/** The choices for a picker param, or `null` for a free-entry param. */
export function paramOptions(
    param: JobParamSpec,
    values: JobFormValues,
    sources: JobsDataSource[],
): JobParamOption[] | null {
    if (param.type === "datasource") {
        return sources.map((s) => ({ value: s.name, label: s.name }));
    }
    if (param.type === "subgroup" || param.type === "column") {
        const source = sources.find((s) => s.name === values[param.options_from ?? ""]);
        if (!source) return [];
        if (param.type === "subgroup") return source.subgroups.map((k) => ({ value: k, label: k }));
        return source.columns.map((c) => ({ value: c.field, label: c.name }));
    }
    return null;
}

function initialValue(param: JobParamSpec, values: JobFormValues, sources: JobsDataSource[]): string {
    const options = paramOptions(param, values, sources);
    if (options === null) return param.default === null ? "" : String(param.default);
    // columns stay unchosen, so the user picks them on purpose
    if (param.type === "column") return "";
    if (options.some((o) => o.value === param.default)) return String(param.default);
    return options[0]?.value ?? "";
}

/** Starting values for a tool. Params are filled in order, so a bound param sees its datasource. */
export function initialValues(tool: JobTool, sources: JobsDataSource[]): JobFormValues {
    const values: JobFormValues = {};
    for (const p of tool.params) {
        values[p.name] = initialValue(p, values, sources);
    }
    return values;
}

/** Set one param; params bound to it through `options_from` go back to their starting value. */
export function setParam(
    tool: JobTool,
    values: JobFormValues,
    name: string,
    value: string,
    sources: JobsDataSource[],
): JobFormValues {
    const next = { ...values, [name]: value };
    for (const p of tool.params) {
        if (p.options_from === name) next[p.name] = initialValue(p, next, sources);
    }
    return next;
}

/** The params to POST. Throws a message for the user when a numeric param does not parse. */
export function toSubmitParams(tool: JobTool, values: JobFormValues): Record<string, JobParamValue> {
    const params: Record<string, JobParamValue> = {};
    for (const p of tool.params) {
        const raw = (values[p.name] ?? "").trim();
        if (p.type === "int" || p.type === "float") {
            // an empty numeric param is left out, so the worker uses its default
            if (raw === "") continue;
            const n = Number(raw);
            if (!Number.isFinite(n)) throw new Error(`${p.label} must be a number`);
            if (p.type === "int" && !Number.isInteger(n)) throw new Error(`${p.label} must be a whole number`);
            params[p.name] = n;
        } else {
            if (raw === "" && PICKERS.includes(p.type)) throw new Error(`Choose a value for ${p.label}`);
            params[p.name] = values[p.name] ?? "";
        }
    }
    return params;
}
