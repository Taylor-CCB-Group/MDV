import { describe, expect, test } from "vitest";
import { RowsAsColsQuery, type RowsAsColslink } from "@/links/link_utils";

function mockLink(): RowsAsColslink {
    return {
        name_column: "gene",
        name: "genes",
        subgroups: {
            gs: { name: "gs", label: "Main", type: "double" },
            layer_a: { name: "layer_a", label: "Layer A", type: "double" },
        },
        observableFields: [],
        initPromise: Promise.resolve(),
    };
}

describe("RowsAsColsQuery", () => {
    test("toJSON always includes canonical subgroupName (default → first key)", () => {
        const link = mockLink();
        const q = new RowsAsColsQuery(link, "genes", 2);
        expect(q.toJSON()).toEqual({
            linkedDsName: "genes",
            maxItems: 2,
            type: "RowsAsColsQuery",
            subgroupName: "gs",
        });
    });

    test("toJSON includes explicit subgroup when set", () => {
        const link = mockLink();
        const q = new RowsAsColsQuery(link, "genes", 4, "layer_a");
        expect(q.toJSON()).toEqual({
            linkedDsName: "genes",
            maxItems: 4,
            type: "RowsAsColsQuery",
            subgroupName: "layer_a",
        });
    });

    test("effectiveSubgroupKey prefers stored subgroup when valid", () => {
        const link = mockLink();
        const q = new RowsAsColsQuery(link, "genes", 1, "layer_a");
        expect(q.effectiveSubgroupKey).toBe("layer_a");
    });

    test("effectiveSubgroupKey falls back to first subgroup when name invalid", () => {
        const link = mockLink();
        const q = new RowsAsColsQuery(link, "genes", 1, "missing");
        expect(q.effectiveSubgroupKey).toBe("gs");
    });

    test("toJSON serializes canonical subgroup when stored key is invalid", () => {
        const link = mockLink();
        const q = new RowsAsColsQuery(link, "genes", 1, "bogus");
        expect(q.toJSON()).toEqual({
            linkedDsName: "genes",
            maxItems: 1,
            type: "RowsAsColsQuery",
            subgroupName: "gs",
        });
    });
});
