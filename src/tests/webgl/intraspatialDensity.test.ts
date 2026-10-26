import { describe, expect, test } from "vitest";
import { contourToDensityProps } from "@/webgl/intraspatialDensity";

type Props = Parameters<typeof contourToDensityProps>[0];

const x = new Float32Array([0, 1, 2, 3]);
const y = new Float32Array([0, 1, 2, 3]);
const indices = [0, 2, 3];

function props(extra: Partial<Props>): Props {
    return {
        id: "c",
        data: indices,
        x,
        y,
        fillOpacity: 0.5,
        contourOpacity: 0.5,
        contourFill: 10000,
        getPosition: (_i: number, { target }: { target: number[] | Float32Array }) => target,
        colorRange: [[255, 0, 0]],
        radiusPixels: 30,
        bandwidth: 1,
        debounce: 1000,
        weightsTextureSize: 256,
        pickable: false,
        ...extra,
    };
}

/** Field-contour style: a fresh closure each call, as each render makes one. */
function fieldProps(field: Float32Array, min: number, max: number) {
    return props({
        getWeight: (i: number) => (field[i] - min) / (max - min),
        weightKey: [field, min, max],
    });
}

describe("contourToDensityProps weights", () => {
    test("1 at the drawn indices, 0 elsewhere, and NaN weights dropped", () => {
        expect(Array.from(contourToDensityProps(props({})).weights ?? [])).toEqual([1, 0, 1, 1]);
        const nan = props({ getWeight: (i: number) => (i === 2 ? Number.NaN : 0.5) });
        expect(Array.from(contourToDensityProps(nan).weights ?? [])).toEqual([0.5, 0, 0, 0.5]);
    });

    test("a keyed weight is reused across new closures, per field, and recomputed when its range changes", () => {
        const a = new Float32Array([0, 1, 2, 4]);
        const b = new Float32Array([4, 2, 1, 0]);
        const wa = contourToDensityProps(fieldProps(a, 0, 4)).weights;
        const wb = contourToDensityProps(fieldProps(b, 0, 4)).weights;
        // fields sharing one index list keep separate entries
        expect(contourToDensityProps(fieldProps(a, 0, 4)).weights).toBe(wa);
        expect(contourToDensityProps(fieldProps(b, 0, 4)).weights).toBe(wb);
        const narrowed = contourToDensityProps(fieldProps(a, 0, 2)).weights;
        expect(narrowed).not.toBe(wa);
        expect(Array.from(narrowed ?? [])).toEqual([0, 0, 1, 2]);
    });
});
