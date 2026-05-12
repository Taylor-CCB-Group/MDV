import { describe, expect, test, vi } from "vitest";

vi.mock("@/react/components/TableChartReactComponent", () => ({
    default: () => null,
}));
vi.mock("@/react/components/BaseReactChart", () => ({
    BaseReactChart: class {},
}));
vi.mock("@/charts/BaseChart", () => ({
    default: class {
        static types: Record<string, unknown> = {};
    },
}));

import { orderSerializedTableParams } from "@/react/components/TableChartReactWrapper";

describe("orderSerializedTableParams", () => {
    test("retains an active-link param when it currently resolves to zero visible fields", () => {
        const serializedQuery = {
            linkedDsName: "genes",
            maxItems: 10,
            type: "RowsAsColsQuery" as const,
        };
        const runtimeQuery = {
            columns: [],
            fields: [] as string[],
            initialize: async () => undefined,
        };

        const result = orderSerializedTableParams(
            [runtimeQuery, "age"],
            [serializedQuery, "age"],
            ["age"],
        );

        expect(result).toEqual([serializedQuery, "age"]);
    });

    test("keeps an inactive active-link entry in place while reordering visible plain columns", () => {
        const serializedQuery = {
            linkedDsName: "genes",
            maxItems: 10,
            type: "RowsAsColsQuery" as const,
        };
        const runtimeQuery = {
            columns: [],
            fields: [] as string[],
            initialize: async () => undefined,
        };

        const result = orderSerializedTableParams(
            ["age", runtimeQuery, "name"],
            ["age", serializedQuery, "name"],
            ["name", "age"],
        );

        expect(result).toEqual(["name", serializedQuery, "age"]);
    });
});
