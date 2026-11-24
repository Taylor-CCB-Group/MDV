import type { Viewport } from "@deck.gl/core";
import { describe, expect, test } from "vitest";
import { textLabelAt } from "@/lib/textLayerHit";

/**
 * A stand-in for deck's TextLayer: two labels, each 200 x 40 font-atlas units centred on its
 * position (as `getBoundingRect` reports a centred label), drawn at size 12 from a 48-unit atlas, so
 * a quarter of that (50 x 10 px), with 8 x 2 px of padding.
 */
function labels(fontSize: number | null = 48) {
    return {
        props: {
            data: [{ position: [0, 0] }, { position: [100, 0] }],
            getPosition: (d: { position: number[] }) => d.position,
            getSize: 12,
            sizeScale: 1,
            backgroundPadding: [8, 2],
        },
        state: fontSize === null ? {} : { fontAtlasManager: { props: { fontSize } } },
        getBoundingRect: () => [-100, -20, 200, 40],
    };
}

/** Data x, y are pixels, offset by (10, 20). */
const viewport = { project: ([x, y]: number[]) => [x + 10, y + 20] } as unknown as Viewport;

describe("which label a press is on", () => {
    test("inside a label's box, padding included", () => {
        expect(textLabelAt(labels(), [10, 20], viewport)).toBe(0);
        expect(textLabelAt(labels(), [10 - 25 - 7, 20 + 5 + 1], viewport)).toBe(0);
        expect(textLabelAt(labels(), [110 + 30, 20], viewport)).toBe(1);
    });

    test("outside every box", () => {
        expect(textLabelAt(labels(), [10 - 25 - 9, 20], viewport)).toBe(-1);
        expect(textLabelAt(labels(), [10, 20 + 5 + 3], viewport)).toBe(-1);
        expect(textLabelAt(labels(), [60, 20], viewport)).toBe(-1);
    });

    test("can't tell before the font atlas is laid out", () => {
        expect(textLabelAt(labels(null), [10, 20], viewport)).toBeNull();
    });
});
