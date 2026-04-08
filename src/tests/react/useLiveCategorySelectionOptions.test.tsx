import { useLiveCategorySelectionValues } from "@/react/components/categorySelectionHooks";
import { act, renderHook } from "@testing-library/react";
import { describe, expect, test, vi } from "vitest";

const mockUseDataStore = vi.fn();

vi.mock("@/react/context", () => ({
    useChart: vi.fn(),
    useDataStore: () => mockUseDataStore(),
}));

describe("useLiveCategorySelectionValues", () => {
    test("refreshes multitext options when the source column gains new categories", () => {
        const listeners: Record<string, (type: string, data?: unknown) => void> = {};
        const tagColumn = {
            datatype: "multitext" as const,
            delimiter: ",",
            field: "__tags",
            stringLength: 1,
            values: ["a"],
        };
        const dataStore = {
            addListener: vi.fn((id: string, listener: (type: string, data?: unknown) => void) => {
                listeners[id] = listener;
            }),
            columnIndex: {
                __tags: tagColumn,
            },
            getColumnValues: vi.fn(),
            removeListener: vi.fn((id: string) => {
                delete listeners[id];
            }),
        };
        mockUseDataStore.mockReturnValue(dataStore);

        const { result } = renderHook(() =>
            useLiveCategorySelectionValues(tagColumn.field),
        );

        expect(result.current).toEqual(["a"]);

        act(() => {
            tagColumn.values.push("a, b", "b");
            Object.values(listeners).forEach((listener) =>
                listener("data_changed", { columns: ["__tags"] }),
            );
        });

        expect(result.current).toEqual(["a", "b"]);
    });
});
