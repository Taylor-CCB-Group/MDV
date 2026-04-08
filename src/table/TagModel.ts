import { loadColumn } from "@/dataloaders/DataLoaderUtil";
import { isColumnLoaded, isColumnOfType } from "@/lib/columnTypeHelpers";
import {
    MULTITEXT_EMPTY_INDEX,
    getMultitextCapacity,
    getMultitextDelimiter,
    getMultitextItemValues,
    getRowMultitextItems,
    inferMultitextEmptyItem,
    isPackedMultitextColumn,
    joinMultitextItems,
    normalizeMultitextItems,
    splitMultitextItems,
} from "@/lib/multitext";
import type { LoadedDataColumn } from "../charts/charts";
import type DataStore from "../datastore/DataStore";
import { DataModel } from "./DataModel";

const COL_TYPE = "multitext" as const;
const MAX_MULTITEXT_VALUES = 65535;

export type TagColumn = LoadedDataColumn<"multitext">;
export type TagSelectionScope = "filtered" | "highlighted";

export type TagAnnotationViewState = {
    tagList: Set<string>;
    tagsInSelection: Set<string>;
    selectionCount: number;
};

export type TagAnnotationWarningFlags = {
    showNoSelectionWarning: boolean;
    showFilteredScopeCoversWholeTableWarning: boolean;
};

type DataChangedEvent = {
    columns?: string[];
};

function getEmptyItem(column: TagColumn) {
    return inferMultitextEmptyItem(column);
}

function getSemanticRowItems(column: TagColumn, rowIndex: number) {
    return getRowMultitextItems(column, rowIndex, {
        emptyItem: getEmptyItem(column),
    });
}

function getSemanticColumnItems(column: TagColumn) {
    return getMultitextItemValues(column, {
        emptyItem: getEmptyItem(column),
    });
}

function getQueryTags(column: TagColumn, tag: string) {
    return splitMultitextItems(tag, getMultitextDelimiter(column))
        .map((item) => item.trim())
        .filter(Boolean)
        .filter((item) => item !== getEmptyItem(column));
}

function getOrAddValueIndex(value: string, column: TagColumn) {
    let index = column.values.indexOf(value);
    if (index !== -1) {
        return index;
    }

    if (column.values.length >= MAX_MULTITEXT_VALUES) {
        throw new Error(
            `multitext column '${column.name}' exceeded ${MAX_MULTITEXT_VALUES} values when adding '${value}'`,
        );
    }

    column.values.push(value);
    index = column.values.length - 1;
    return index;
}

type RowWritePlan =
    | { kind: "packed"; itemsToWrite: string[] }
    | { kind: "legacy"; encodedValue: string };

function computeRowWritePlan(column: TagColumn, items: string[]): RowWritePlan {
    const normalizedItems = normalizeMultitextItems(items);
    const capacity = getMultitextCapacity(column);

    if (isPackedMultitextColumn(column)) {
        const emptyItem = getEmptyItem(column);
        const itemsToWrite =
            normalizedItems.length === 0 && emptyItem
                ? [emptyItem]
                : normalizedItems;

        if (itemsToWrite.length > capacity) {
            throw new Error(
                `multitext column '${column.name}' can only store ${capacity} values per row`,
            );
        }

        return { kind: "packed", itemsToWrite };
    }

    const encodedValue =
        normalizedItems.length === 0
            ? ""
            : joinMultitextItems(normalizedItems, getMultitextDelimiter(column));

    return { kind: "legacy", encodedValue };
}

function valueStringsFromRowWritePlan(plan: RowWritePlan): string[] {
    if (plan.kind === "packed") {
        return plan.itemsToWrite;
    }
    return plan.encodedValue === "" ? [] : [plan.encodedValue];
}

/**
 * Ensures the column dictionary can accept every distinct new string in one shot.
 * Call before mutating row buffers or column.values so failures cannot leave partial writes.
 */
function assertCanAddNewMultitextValues(
    column: TagColumn,
    valueStrings: readonly string[],
): void {
    const distinctNew = new Set<string>();
    for (const v of valueStrings) {
        if (column.values.indexOf(v) === -1) {
            distinctNew.add(v);
        }
    }
    if (column.values.length + distinctNew.size > MAX_MULTITEXT_VALUES) {
        throw new Error(
            `multitext column '${column.name}' exceeded ${MAX_MULTITEXT_VALUES} values`,
        );
    }
}

function collectPreflightValueStringsForUpdate(
    tagToChange: string,
    column: TagColumn,
    rowIndices: Iterable<number>,
    tagValue: boolean,
    rowState?: Map<number, string[]>,
): string[] {
    const tags = splitMultitextItems(
        tagToChange,
        getMultitextDelimiter(column),
    );

    if (tags.length > 1) {
        const acc: string[] = [];
        const sim = rowState ?? new Map<number, string[]>();
        for (const tag of tags) {
            acc.push(
                ...collectPreflightValueStringsForUpdate(
                    tag,
                    column,
                    rowIndices,
                    tagValue,
                    sim,
                ),
            );
        }
        return acc;
    }

    const tag = tags[0]?.trim() || tagToChange.trim();
    if (!tag || tag === getEmptyItem(column)) {
        return [];
    }

    const pending: string[] = [];
    const sim = rowState ?? new Map<number, string[]>();
    for (const rowIndex of rowIndices) {
        const currentItems =
            sim.get(rowIndex) ?? getSemanticRowItems(column, rowIndex);
        const hasTag = currentItems.includes(tag);
        if (tagValue ? hasTag : !hasTag) {
            continue;
        }

        const nextItems = tagValue
            ? [...currentItems, tag]
            : currentItems.filter((item) => item !== tag);
        sim.set(rowIndex, nextItems);
        pending.push(
            ...valueStringsFromRowWritePlan(
                computeRowWritePlan(column, nextItems),
            ),
        );
    }
    return pending;
}

function clearRow(column: TagColumn, rowIndex: number) {
    const capacity = getMultitextCapacity(column);
    const baseIndex = rowIndex * capacity;
    for (let offset = 0; offset < capacity; offset += 1) {
        column.data[baseIndex + offset] = MULTITEXT_EMPTY_INDEX;
    }
}

function writeRowItems(column: TagColumn, rowIndex: number, items: string[]) {
    const plan = computeRowWritePlan(column, items);
    assertCanAddNewMultitextValues(column, valueStringsFromRowWritePlan(plan));

    const capacity = getMultitextCapacity(column);
    const baseIndex = rowIndex * capacity;

    clearRow(column, rowIndex);
    if (plan.kind === "packed") {
        plan.itemsToWrite.forEach((item, offset) => {
            column.data[baseIndex + offset] = getOrAddValueIndex(item, column);
        });
        return;
    }
    if (plan.encodedValue !== "") {
        column.data[baseIndex] = getOrAddValueIndex(plan.encodedValue, column);
    }
}

function updateSelectedRows(
    tagToChange: string,
    column: TagColumn,
    rowIndices: Iterable<number>,
    dataStore: DataStore,
    notify = true,
    tagValue = true,
): boolean {
    const tags = splitMultitextItems(
        tagToChange,
        getMultitextDelimiter(column),
    );

    assertCanAddNewMultitextValues(
        column,
        collectPreflightValueStringsForUpdate(
            tagToChange,
            column,
            rowIndices,
            tagValue,
        ),
    );

    if (tags.length > 1) {
        const changed = tags
            .map((tag) =>
                updateSelectedRows(tag, column, rowIndices, dataStore, false, tagValue),
            )
            .some(Boolean);

        if (notify && changed) {
            dataStore.dataChanged([column.name]);
        }
        return changed;
    }

    const tag = tags[0]?.trim() || tagToChange.trim();
    if (!tag || tag === getEmptyItem(column)) {
        return false;
    }

    let changed = false;
    for (const rowIndex of rowIndices) {
        const currentItems = getSemanticRowItems(column, rowIndex);
        const hasTag = currentItems.includes(tag);
        if (tagValue ? hasTag : !hasTag) {
            continue;
        }

        const nextItems = tagValue
            ? [...currentItems, tag]
            : currentItems.filter((item) => item !== tag);
        writeRowItems(column, rowIndex, nextItems);
        changed = true;
    }

    if (notify && changed) {
        dataStore.dataChanged([column.name]);
    }
    return changed;
}

/**
 * TagModel manages item-level annotations stored in a multitext column.
 *
 * It understands both legacy single-slot columns whose stored values are joined
 * strings such as `"a; b"` and packed multitext columns such as the gates column.
 */
export default class TagModel {
    tagColumn: TagColumn;
    readonly dataStore: DataStore;
    readonly dataModel: DataModel;
    isReady = false;
    private listeners: (() => void)[] = [];
    private readonly selectionScope: TagSelectionScope;
    private readonly listenerId: string;

    private constructor(
        dataStore: DataStore,
        columnName = "__tags",
        selectionScope: TagSelectionScope = "filtered",
    ) {
        this.dataStore = dataStore;
        this.selectionScope = selectionScope;
        this.listenerId = `tag-model-${columnName}-${Math.random().toString(36).slice(2)}`;

        const column = dataStore.columnIndex[columnName];
        if (!column) {
            const columnSpec = {
                name: columnName,
                field: columnName,
                datatype: COL_TYPE,
                editable: true,
                delimiter: ";" as const,
                stringLength: 1,
                values: [] as string[],
            };
            const buffer = new SharedArrayBuffer(
                dataStore.size * Uint16Array.BYTES_PER_ELEMENT,
            );
            const data = new Uint16Array(buffer);
            data.fill(MULTITEXT_EMPTY_INDEX);
            dataStore.addColumn(columnSpec, buffer);
            this.tagColumn = dataStore.columnIndex[columnName] as TagColumn;
            this.isReady = true;
        } else {
            if (isColumnOfType(column, COL_TYPE) && column.editable === false) {
                throw new Error(`column '${columnName}' is not editable`);
            }
            this.tagColumn = column as TagColumn;
        }

        this.dataModel = new DataModel(dataStore, { autoupdate: true });
        this.dataModel.updateModel();
        this.dataModel.addListener("tag", () => {
            if (this.selectionScope === "filtered") {
                this.callListeners();
            }
        });
        this.dataStore.addListener(this.listenerId, (type, data?: DataChangedEvent) => {
            if (type === "data_highlighted" && this.selectionScope === "highlighted") {
                this.callListeners();
                return;
            }

            if (
                type === "data_changed"
                && Array.isArray(data?.columns)
                && data.columns.includes(this.tagColumn.name)
            ) {
                if (this.selectionScope === "filtered") {
                    this.dataModel.updateModel();
                }
                this.callListeners();
            }
        });
    }

    static async create(
        dataStore: DataStore,
        columnName = "__tags",
        selectionScope: TagSelectionScope = "filtered",
    ): Promise<TagModel> {
        const model = new TagModel(dataStore, columnName, selectionScope);

        if (!model.isReady) {
        try {
            const column = await loadColumn(dataStore.name, columnName);
            if (!isColumnOfType(column, COL_TYPE)) {
                throw new Error(
                    `column '${columnName}' is not of type '${COL_TYPE}'`,
                );
            }
            if (!isColumnLoaded(column)) {
                throw new Error("expected column to be loaded");
            }
            if (column.editable === false) {
                throw new Error(`column '${columnName}' is not editable`);
            }
            model.tagColumn = column;
            model.isReady = true;
            model.callListeners();
        } catch (error) {
            model.dispose();
            throw error;
            }
        }

        return model;
    }

    private getSelectedRowIndices() {
        if (this.selectionScope === "highlighted") {
            return this.dataStore.getHighlightedData()?.slice() || [];
        }
        return this.dataModel.data;
    }

    getSelectionLength() {
        return this.getSelectedRowIndices().length;
    }

    getAnnotationViewState(): TagAnnotationViewState {
        return {
            tagList: this.getTags(),
            tagsInSelection: this.getTagsInSelection(),
            selectionCount: this.getSelectionLength(),
        };
    }

    /**
     * Warning flags for annotation UI: empty selection vs filtered scope covering the full table
     * (no narrowing filter, or selection length equals datastore row count).
     */
    getAnnotationWarningFlags(): TagAnnotationWarningFlags {
        const selectionCount = this.getSelectionLength();
        return {
            showNoSelectionWarning: selectionCount === 0,
            showFilteredScopeCoversWholeTableWarning:
                selectionCount > 0 && this.isFilteredScopeAnnotatingEntireTable(),
        };
    }

    private isFilteredScopeAnnotatingEntireTable(): boolean {
        if (this.selectionScope !== "filtered") {
            return false;
        }

        if (typeof this.dataStore.isFiltered === "function") {
            return !this.dataStore.isFiltered();
        }

        const totalRows = this.dataStore.size ?? this.getSelectionLength();
        return this.getSelectionLength() === totalRows;
    }

    private callListeners() {
        this.listeners.forEach((callback) => callback());
    }

    addListener(callback: () => void) {
        this.listeners.push(callback);
        return callback;
    }

    removeListener(callback: () => void) {
        const index = this.listeners.indexOf(callback);
        if (index !== -1) {
            this.listeners.splice(index, 1);
        }
    }

    dispose() {
        this.dataStore.removeListener(this.listenerId);
        this.dataModel.dispose();
        this.listeners = [];
    }

    setTag(tag: string, tagValue = true) {
        updateSelectedRows(
            tag,
            this.tagColumn,
            this.getSelectedRowIndices(),
            this.dataStore,
            true,
            tagValue,
        );
        this.dataModel.updateModel();
        this.callListeners();
    }

    getTags(): Set<string> {
        return new Set(getSemanticColumnItems(this.tagColumn));
    }

    entireSelectionHasTag(tag: string): boolean {
        const tags = splitMultitextItems(tag, getMultitextDelimiter(this.tagColumn));
        if (tags.length > 1) {
            return tags.every((item) => this.entireSelectionHasTag(item));
        }

        const needle = tags[0]?.trim() || tag.trim();
        if (!needle) {
            return false;
        }

        return !Array.from(this.getSelectedRowIndices()).some(
            (rowIndex) => !getSemanticRowItems(this.tagColumn, rowIndex).includes(needle),
        );
    }

    getTagsInSelection(): Set<string> {
        const tags = new Set<string>();
        for (const rowIndex of this.getSelectedRowIndices()) {
            getSemanticRowItems(this.tagColumn, rowIndex).forEach((tag) =>
                tags.add(tag),
            );
        }
        return tags;
    }

    getMatchingRowIndices(tag: string) {
        const queryTags = getQueryTags(this.tagColumn, tag);
        if (queryTags.length === 0) {
            return [];
        }

        return Array.from(this.getSelectedRowIndices()).filter((rowIndex) => {
            const rowItems = getSemanticRowItems(this.tagColumn, rowIndex);
            return queryTags.every((queryTag) => rowItems.includes(queryTag));
        });
    }
}
