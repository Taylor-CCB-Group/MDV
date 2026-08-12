import type { RenderStack } from "@spatialdata/layers";
import { autorun, observable, runInAction } from "mobx";
import { describe, expect, test } from "vitest";

import {
    mdvFieldSpecsOf,
    projectMdvFieldSpecs,
    withoutMdvFieldSpecs,
} from "@/react/spatialdata/field_spec_projection";

/**
 * An "active link" is a column choice that keeps moving: its columns are whatever
 * the linked datasource currently has selected. The viewer has no vocabulary for
 * that — `fillColorByColumn.columnName` is a plain string — so MDV keeps the query
 * on its own side and writes the concrete name from it. These tests are about that
 * projection, which is the only thing standing between a live link and a layer that
 * silently never changes.
 */

/** Stands in for a `RowsAsColsQuery`: `flattenFields` reads `.fields` off anything non-string. */
function stubQuery(fields: string[]) {
    return observable({ fields }, {}, { deep: false });
}

function stackWith(props: Record<string, unknown>): RenderStack {
    return observable({
        entries: [
            {
                kind: "spatial",
                id: "spatialdata-shapes-cells-#1",
                visible: true,
                source: { elementType: "shapes", elementKey: "cells" },
                props,
            },
        ],
    }) as unknown as RenderStack;
}

describe("projecting MDV's field specs onto the viewer's fields", () => {
    test("writes the column a query currently resolves to", () => {
        const stack = stackWith({
            fillColorByColumn: { columnName: "cell_type", mode: "auto" },
            mdvFieldSpecs: { fillColorByColumn: stubQuery(["gs|CD163 (gs)|7"]) },
        });

        expect(projectMdvFieldSpecs(stack)).toBe(true);
        expect(stack.entries[0].props.fillColorByColumn).toEqual({
            columnName: "gs|CD163 (gs)|7",
            mode: "auto",
        });
    });

    test("follows the link when it resolves somewhere else", () => {
        const query = stubQuery(["gs|CD163 (gs)|7"]);
        const stack = stackWith({ mdvFieldSpecs: { fillColorByColumn: query } });

        const seen: unknown[] = [];
        const dispose = autorun(() => {
            projectMdvFieldSpecs(stack);
            seen.push((stack.entries[0].props as Record<string, unknown>).fillColorByColumn);
        });

        // What selecting a different gene in the linked table amounts to.
        runInAction(() => {
            query.fields = ["gs|CD31 (gs)|11"];
        });
        dispose();

        expect(seen.at(-1)).toEqual({ columnName: "gs|CD31 (gs)|11" });
    });

    test("reports no change when the resolution is the same, so nothing downstream is nudged", () => {
        const stack = stackWith({
            fillColorByColumn: { columnName: "cell_type", mode: "auto" },
            mdvFieldSpecs: { fillColorByColumn: "cell_type" },
        });

        expect(projectMdvFieldSpecs(stack)).toBe(false);
    });

    test("leaves the previous column alone while the link is still initialising", () => {
        // A query that resolves to nothing has not been answered yet — it is not the
        // user asking for no colour. Clearing here would flash the canvas through
        // uncoloured on every reload.
        const stack = stackWith({
            fillColorByColumn: { columnName: "cell_type", mode: "auto" },
            mdvFieldSpecs: { fillColorByColumn: stubQuery([]) },
        });

        expect(projectMdvFieldSpecs(stack)).toBe(false);
        expect(stack.entries[0].props.fillColorByColumn).toEqual({
            columnName: "cell_type",
            mode: "auto",
        });
    });

    test("projects every field of a multi-column tooltip spec", () => {
        const stack = stackWith({
            tooltipFields: ["cell_type"],
            mdvFieldSpecs: { tooltipFields: ["cell_type", stubQuery(["gs|CD163 (gs)|7"])] },
        });

        expect(projectMdvFieldSpecs(stack)).toBe(true);
        expect(stack.entries[0].props.tooltipFields).toEqual(["cell_type", "gs|CD163 (gs)|7"]);
    });

    test("leaves a layer with no specs to its own devices", () => {
        // Every layer saved before this existed, and every one whose column was set
        // some other way: the concrete fields stand on their own.
        const stack = stackWith({ fillColorByColumn: { columnName: "cell_type", mode: "auto" } });

        expect(projectMdvFieldSpecs(stack)).toBe(false);
        expect(stack.entries[0].props.fillColorByColumn).toEqual({
            columnName: "cell_type",
            mode: "auto",
        });
    });
});

describe("the viewer boundary", () => {
    test("takes the specs off, because a live query is not layer config", () => {
        const config = { type: "shapes", elementKey: "cells", mdvFieldSpecs: { fillColorByColumn: "x" } };

        expect(withoutMdvFieldSpecs(config)).toEqual({ type: "shapes", elementKey: "cells" });
    });

    test("passes a config that never had them through untouched", () => {
        const config = { type: "shapes", elementKey: "cells" };

        expect(withoutMdvFieldSpecs(config)).toBe(config);
    });

    test("reads nothing out of a value that is not a specs object", () => {
        expect(mdvFieldSpecsOf(undefined)).toBeUndefined();
        expect(mdvFieldSpecsOf({ mdvFieldSpecs: "not an object" })).toBeUndefined();
        expect(mdvFieldSpecsOf({ mdvFieldSpecs: ["not an object either"] })).toBeUndefined();
    });
});
