import { autorun } from "mobx";
import { afterEach, beforeEach, describe, expect, test, vi } from "vitest";
import DataStore from "@/datastore/DataStore";
import { getFieldName, RowsAsColsQuery } from "@/links/link_utils";
import { type InProcessFilteredIndexWorker, stubWorker } from "./dimensionHarness";

/**
 * Rows-as-columns links against real DataStores: two dataSources (e.g. cells and a grid of spots) that each link
 * to the same genes dataSource, using the same link name.
 */

const GENES = ["A", "B", "C", "D"];
let genes: DataStore;
let cells: DataStore;
let grid: DataStore;
let workers: InProcessFilteredIndexWorker[];

/** a dataSource with a rows-as-columns link to `genes` - the link config is not shared between dataSources */
function linkedStore(name: string) {
    return new DataStore(3, {
        name,
        size: 3,
        columns: [],
        links: {
            genes: {
                rows_as_columns: {
                    name: "Gene Expr",
                    name_column: "name",
                    subgroups: { gs: { name: "gs", label: "Gene scores", type: "sparse" } },
                },
            },
        },
    } as never);
}

const query = (ds: DataStore) =>
    RowsAsColsQuery.fromSerialized(ds, { linkedDsName: "genes", maxItems: 10, type: "RowsAsColsQuery" });

/** Lets promises settle and has the filtered index workers answer, until there is nothing left to do. */
async function settle() {
    for (let quiet = 0; quiet < 3; ) {
        await new Promise((resolve) => setTimeout(resolve, 0));
        const waiting = workers.filter((w) => w.queue.length > 0);
        for (const w of waiting) w.runNext();
        quiet = waiting.length ? 0 : quiet + 1;
    }
}

const fieldsFor = (...rows: number[]) => rows.map((i) => getFieldName("gs", GENES[i], i));

beforeEach(() => {
    workers = stubWorker({ filteredIndex: true }).filteredIndexWorkers;
    vi.spyOn(console, "log").mockImplementation(() => {});
    genes = new DataStore(GENES.length, {
        name: "genes",
        size: GENES.length,
        columns: [{ field: "name", name: "name", datatype: "text" }],
    } as never);
    genes.setColumnData("name", GENES);
    cells = linkedStore("cells");
    grid = linkedStore("grid");
    vi.stubGlobal("mdv", {
        chartManager: {
            dataSources: [genes, cells, grid].map((dataStore) => ({ name: dataStore.name, dataStore })),
            loadColumnSetAsync: async () => {},
        },
    });
});
afterEach(() => {
    vi.restoreAllMocks();
    vi.unstubAllGlobals();
});

describe("two dataSources linked to the same rows-as-columns dataSource", () => {
    test("both links follow a highlight", async () => {
        const fromCells = query(cells);
        const fromGrid = query(grid);
        await settle();
        expect(fromCells.fields).toEqual(fieldsFor(0, 1, 2, 3));
        expect(fromGrid.fields).toEqual(fieldsFor(0, 1, 2, 3));

        genes.dataHighlighted([2], {});
        await settle();

        expect(fromGrid.fields).toEqual(fieldsFor(2));
        expect(fromCells.fields).toEqual(fieldsFor(2));
    });

    test("both links follow a filter", async () => {
        const fromCells = query(cells);
        const fromGrid = query(grid);
        await settle();

        genes.getDimension("category_dimension").filter("filterCategories", ["name"], ["B", "D"]);
        await settle();

        expect(fromGrid.fields).toEqual(fieldsFor(1, 3));
        expect(fromCells.fields).toEqual(fieldsFor(1, 3));
    });
});

describe("a link's fields over time", () => {
    /** every value `query.fields` takes from now on, as something observing the link would see them */
    function recordFields(q: RowsAsColsQuery) {
        const history: string[][] = [];
        const stop = autorun(() => {
            history.push(q.fields);
        });
        return { history, stop };
    }

    test("a link that initialises when rows are already highlighted only ever exposes those rows", async () => {
        // as when a view is loaded: its highlight is restored before any chart resolves the link
        genes.dataHighlighted([2], {});
        const { history, stop } = recordFields(query(cells));
        await settle();
        stop();

        expect(history.filter((fields) => fields.length > 0)).toEqual([fieldsFor(2)]);
    });

    test("an event that leaves the rows as they were is not a change to the link's fields", async () => {
        genes.dataHighlighted([2], {});
        const { history, stop } = recordFields(query(cells));
        await settle();
        const before = history.length;

        // the highlight is announced again once a view's charts exist, for those that only show it on the event
        genes.dataHighlighted([2], {});
        genes.triggerFilter();
        await settle();
        stop();

        expect(history.length).toBe(before);
    });
});
