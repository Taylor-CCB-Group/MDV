import { MULTITEXT_EMPTY_INDEX } from "@/lib/multitext";
import { getCategoryFilterIndices } from "@/lib/categoryFilterUtils";
import { describe, expect, test } from "vitest";

describe("getCategoryFilterIndices", () => {
    test("matches multitext rows that include the selected item", () => {
        const contourParameter = {
            name: "Tags",
            field: "__tags",
            datatype: "multitext" as const,
            values: ["a, b", "b", "c"],
            data: new Uint8Array([0, 1, 2]),
            getValue: () => "",
            stringLength: 1,
        };

        expect(
            getCategoryFilterIndices(
                new Uint32Array([0, 1, 2]),
                contourParameter,
                ["b"],
            ),
        ).toEqual(new Uint32Array([0, 1]));
    });

    test("preserves exact raw-value matching for legacy multitext selections", () => {
        const contourParameter = {
            name: "Tags",
            field: "__tags",
            datatype: "multitext" as const,
            values: ["a, b", "b", "c"],
            data: new Uint8Array([0, 1, 2]),
            getValue: () => "",
            stringLength: 1,
        };

        expect(
            getCategoryFilterIndices(
                new Uint32Array([0, 1, 2]),
                contourParameter,
                "a, b",
            ),
        ).toEqual(new Uint32Array([0]));
    });

    test("matches packed multitext rows containing the selected item", () => {
        const contourParameter = {
            name: "Tags",
            field: "__tags",
            datatype: "multitext" as const,
            values: ["N/A", "a", "b", "c"],
            data: new Uint16Array([
                1, 2,
                2, MULTITEXT_EMPTY_INDEX,
                3, MULTITEXT_EMPTY_INDEX,
            ]),
            getValue: () => "",
            stringLength: 2,
            delimiter: ";",
        };

        expect(
            getCategoryFilterIndices(
                new Uint32Array([0, 1, 2]),
                contourParameter,
                ["b"],
            ),
        ).toEqual(new Uint32Array([0, 1]));
    });
});
