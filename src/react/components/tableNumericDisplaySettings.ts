import { action, observable, untracked } from "mobx";
import type { TableNumericDisplay } from "@/charts/schemas/ChartConfigSchema";
import { type FieldSpecs, flattenFields } from "@/lib/columnTypeHelpers";
import { g } from "@/lib/utils";
import type { AnyGuiSpec, DataColumn, DataType, GuiSpec } from "../../charts/charts";
import { type ColumnDisplayConfig, isNumericDisplayColumn, withColumnDisplay } from "./tableNumericDisplay";

type TableDisplayConfig = { param: FieldSpecs; column_display?: ColumnDisplayConfig };
type ColumnIndex = Record<string, DataColumn<DataType> | undefined>;

const MODES: [string, TableNumericDisplay["mode"]][] = [
    ["Number", "number"],
    ["Bar", "bar"],
    ["Diverging", "diverging"],
];

function isMode(value: string): value is TableNumericDisplay["mode"] {
    return MODES.some(([, mode]) => mode === value);
}

const getLimit = (display: TableNumericDisplay | undefined) =>
    display?.mode === "diverging" ? display.limit : undefined;
const limitToText = (limit: number | undefined) => (limit === undefined ? "" : String(limit));
const textToLimit = (text: string) => (text.trim() === "" ? undefined : Number(text));
const isValidLimit = (limit: number | undefined) => limit !== undefined && Number.isFinite(limit) && limit > 0;

/**
 * The controls of one column. Their values are getters reading the config, so they follow
 * changes made elsewhere, such as from the column's header menu.
 */
function getColumnSettings(config: TableDisplayConfig, field: string, name: string) {
    const getSaved = () => config.column_display?.[field];
    // what is in the range box, and the display it was typed against
    const draft = observable.box(
        untracked(() => ({ text: limitToText(getLimit(getSaved())), display: getSaved() })),
        { deep: false },
    );
    const setDraft = (text: string) => draft.set({ text, display: getSaved() });

    const getRangeText = () => {
        const saved = getSaved();
        const { text, display } = draft.get();
        // a range saved since the box was typed in replaces it, unless it already means the same
        if (saved === display || saved?.mode !== "diverging" || textToLimit(text) === saved.limit) return text;
        return limitToText(saved.limit);
    };

    const write = action((mode: TableNumericDisplay["mode"], limit?: number) => {
        const display: TableNumericDisplay =
            mode === "diverging" && limit !== undefined ? { mode, limit } : { mode };
        // replaced rather than mutated, so the table sees a new object
        config.column_display = withColumnDisplay(config.column_display, field, display);
    });

    const folder = g({
        type: "folder",
        label: name,
        current_value: [
            g({
                type: "radiobuttons",
                label: "Display",
                get current_value() {
                    return getSaved()?.mode ?? "number";
                },
                // the value is whatever `func` writes to the config
                set current_value(_value: string) {},
                choices: MODES,
                func: action((value: string) => {
                    if (!isMode(value)) return;
                    // remember the range, so it survives switching away from "diverging"
                    if (getSaved()?.mode === "diverging") setDraft(getRangeText());
                    const limit = textToLimit(getRangeText());
                    write(value, isValidLimit(limit) ? limit : undefined);
                }),
            }),
            g({
                type: "text",
                label: "Diverging range ± (blank: data range)",
                get current_value() {
                    return getRangeText();
                },
                set current_value(text: string) {
                    setDraft(text);
                },
                func: (value: string) => {
                    const limit = textToLimit(value);
                    // ignore partial input such as "0." while it is being typed
                    if (limit !== undefined && !isValidLimit(limit)) return;
                    if (getSaved()?.mode === "diverging") write("diverging", limit);
                },
            }),
        ],
    });
    // names need not be unique, so identify the folder by its field
    return Object.assign(folder, { _stableId: `column-display:${field}` });
}

/**
 * A settings folder with the display options of each numeric column in the table.
 *
 * It follows the table through getters rather than reactions with `_disposers`: those are
 * disposed as soon as the dialog mounts in React's strict mode.
 */
export function getNumericDisplaySettings(config: TableDisplayConfig, columnIndex: ColumnIndex): GuiSpec<"folder"> {
    return observable(
        g({
            type: "folder",
            label: "Numeric Column Display",
            // computed from the table's columns, which can change while the dialog is open
            get current_value(): AnyGuiSpec[] {
                return flattenFields(config.param).flatMap((field) => {
                    const column = columnIndex[field];
                    return column && isNumericDisplayColumn(column)
                        ? [getColumnSettings(config, field, column.name)]
                        : [];
                });
            },
        }),
    );
}
