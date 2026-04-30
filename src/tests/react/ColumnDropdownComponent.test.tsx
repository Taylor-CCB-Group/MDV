import { getSelectableColumns } from "@/react/components/columnDropdownUtils";
import { describe, expect, test } from "vitest";

describe("getSelectableColumns", () => {
    test("includes multitext columns for contour-style categorical pickers", () => {
        const columns = [
            {
                datatype: "text",
                field: "cell_type",
                name: "Cell Type",
            },
            {
                datatype: "multitext",
                field: "__tags",
                name: "__tags",
            },
            {
                datatype: "double",
                field: "score",
                name: "Score",
            },
        ] as any;

        expect(
            getSelectableColumns(columns, ["text", "text16", "multitext"]).map(
                (column) => column.field,
            ),
        ).toEqual(["cell_type", "__tags"]);
    });

    test("re-evaluates against the current columns array contents", () => {
        const columns = [
            {
                datatype: "text",
                field: "cell_type",
                name: "Cell Type",
            },
        ] as any[];

        expect(
            getSelectableColumns(columns, ["text", "text16", "multitext"]).map(
                (column) => column.field,
            ),
        ).toEqual(["cell_type"]);

        columns.push({
            datatype: "multitext",
            field: "__tags",
            name: "__tags",
        });

        expect(
            getSelectableColumns(columns, ["text", "text16", "multitext"]).map(
                (column) => column.field,
            ),
        ).toEqual(["cell_type", "__tags"]);
    });
});
