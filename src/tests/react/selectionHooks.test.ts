import {
    getHighlightModifierState,
    getNextHighlightedRows,
} from "@/react/selectionHooks";
import { describe, expect, test } from "vitest";

describe("selectionHooks highlight helpers", () => {
    test("replaces highlights by default", () => {
        expect(getNextHighlightedRows([3, 4], [1, 2])).toEqual([3, 4]);
    });

    test("adds highlights when shift is pressed", () => {
        expect(
            getNextHighlightedRows([2, 3], [1, 2], { shiftKey: true }),
        ).toEqual([1, 2, 3]);
    });

    test("toggles highlights when ctrl or meta is pressed", () => {
        expect(
            getNextHighlightedRows([2, 3], [1, 2], { ctrlKey: true }),
        ).toEqual([1, 2, 3]);
        expect(
            getNextHighlightedRows([2, 3], [1, 2], { metaKey: true }),
        ).toEqual([1, 2, 3]);
    });

    test("removes the whole target group when all toggled rows are already highlighted", () => {
        expect(
            getNextHighlightedRows([2, 3], [1, 2, 3], { ctrlKey: true }),
        ).toEqual([1]);
    });

    test("gives toggle precedence over add when both modifier families are pressed", () => {
        expect(
            getHighlightModifierState({ ctrlKey: true, shiftKey: true }),
        ).toEqual({
            add: false,
            toggle: true,
        });
        expect(
            getNextHighlightedRows([2, 3], [1, 2], {
                ctrlKey: true,
                shiftKey: true,
            }),
        ).toEqual([1, 2, 3]);
    });
});
