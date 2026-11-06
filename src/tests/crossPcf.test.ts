import type { SimpleFeatureCollection } from "@deck.gl-community/editable-layers";
import { Matrix4 } from "@math.gl/core";
import { domainContains } from "intraspatial/spatial/domain";
import { describe, expect, test } from "vitest";
import type { Gate } from "@/react/gates/types";
import {
    type AnalysisFrame,
    frameFromMatrix,
    identityFrame,
    regionsScaleFrame,
    withPhysicalUnit,
} from "@/spatial_stats/analysisFrame";
import { crossPcfInDomain, gatherCells } from "@/spatial_stats/crossPcfCompute";
import {
    featureRings,
    gateDomainSource,
    polygonGeometryDomain,
    selectionDomainSource,
} from "@/spatial_stats/domainSource";
import { regionImageKey, type SpatialDataTransforms, spatialDataRegionFrame } from "@/spatial_stats/spatialDataFrame";

/** A 400×400 square with a 100×100 hole, in data-column units (pixels). */
const holed: SimpleFeatureCollection = {
    type: "FeatureCollection",
    features: [
        {
            type: "Feature",
            properties: {},
            geometry: {
                type: "Polygon",
                coordinates: [
                    [
                        [0, 0],
                        [400, 0],
                        [400, 400],
                        [0, 400],
                        [0, 0],
                    ],
                    [
                        [150, 150],
                        [250, 150],
                        [250, 250],
                        [150, 250],
                        [150, 150],
                    ],
                ],
            },
        },
    ],
};

/** Anisotropic pixels plus an offset: what a real NGFF scale+translation sequence gives. */
const aniso: AnalysisFrame = frameFromMatrix(new Matrix4().translate([10, -5, 0]).scale([0.5, 0.25, 1]), "µm", "test");

/** Stand-in for sd.js: `images[key].getTransformation(cs)` → `{ toMatrix, output.axes }`. */
function fakeStore(matrix: Matrix4, unit?: string): SpatialDataTransforms {
    return {
        images: {
            morphology_focus: {
                getTransformation: (cs) =>
                    cs === "global"
                        ? {
                              ok: true,
                              value: {
                                  toMatrix: () => matrix,
                                  output: {
                                      name: "global",
                                      axes: [
                                          { name: "x", unit },
                                          { name: "y", unit },
                                      ],
                                  },
                              },
                          }
                        : { ok: false, error: new Error(`no ${cs}`) },
            },
        },
    };
}

const region = {
    spatial: { coordinate_system: "global", file: "sdata_X" },
    viv_image: { file: "sdata_X/images/morphology_focus" },
};

describe("analysis frame", () => {
    test("region image key comes from the converter's viv_image path", () => {
        expect(regionImageKey(region)).toBe("morphology_focus");
        expect(regionImageKey({ ...region, viv_image: { file: "other.ome.tiff" } })).toBeUndefined();
    });

    test("SpatialData image transform → µm, converting the coordinate system's unit", () => {
        const px = new Matrix4().scale([0.2125, 0.2125, 1]);
        const um = spatialDataRegionFrame(fakeStore(px, "micrometer"), region);
        if ("error" in um) throw new Error(um.error);
        expect(um.unit).toBe("µm");
        expect(um.toAnalysis[0]).toBeCloseTo(0.2125);

        const mm = spatialDataRegionFrame(fakeStore(new Matrix4().scale([0.001, 0.002, 1]), "millimeter"), region);
        if ("error" in mm) throw new Error(mm.error);
        expect(mm.toAnalysis[0]).toBeCloseTo(1);
        expect(mm.toAnalysis[5]).toBeCloseTo(2);
    });

    test("an unstated unit is reported, not assumed to be µm", () => {
        const f = spatialDataRegionFrame(fakeStore(new Matrix4().identity(), "unit"), region);
        if ("error" in f) throw new Error(f.error);
        expect(f.unit).toContain("unstated");
    });

    test("missing coordinate system or image is an error, not a guess", () => {
        const store = fakeStore(new Matrix4().identity(), "micrometer");
        expect(
            spatialDataRegionFrame(store, { ...region, spatial: { ...region.spatial, coordinate_system: "x" } }),
        ).toHaveProperty("error");
        expect(spatialDataRegionFrame({}, region)).toHaveProperty("error");
    });

    test("regions.scale stopgap stays available and labelled", () => {
        expect(regionsScaleFrame({ scale: 2, scale_unit: "µm" }).provenance).toBe("regions.scale");
        expect(regionsScaleFrame(undefined).provenance).toBe("identity");
    });
});

describe("domain sources", () => {
    test("keeps every ring and maps it through the frame", () => {
        expect(featureRings(holed)).toHaveLength(2);
        const d = polygonGeometryDomain(holed, aniso);
        expect(d.bbox).toEqual([10, -5, 210, 95]);
        // Centre is in the hole; a point between the hole and the outer edge is inside.
        expect(domainContains(d, 10 + 200 * 0.5, -5 + 200 * 0.25)).toBe(false);
        expect(domainContains(d, 10 + 50 * 0.5, -5 + 50 * 0.25)).toBe(true);
    });

    test("a live selection and a saved gate with the same geometry give the same domain", () => {
        const gate: Gate = {
            id: "g",
            name: "g",
            geometry: holed,
            columns: ["x", "y"],
            createdAt: 0,
            labelPosition: [0, 0],
            region: "r1",
        };
        const sel = selectionDomainSource({ selectionFeatureCollection: holed, param: ["x", "y"], region: "r1" });
        const fromGate = gateDomainSource(gate);
        expect(sel?.columns).toEqual(fromGate.columns);
        expect(sel?.region).toBe("r1");
        expect(sel?.toDomain(aniso)).toEqual(fromGate.toDomain(aniso));
    });

    test("an empty selection is no source", () => {
        expect(
            selectionDomainSource({
                selectionFeatureCollection: { type: "FeatureCollection", features: [] },
                param: ["x", "y"],
            }),
        ).toBeUndefined();
    });
});

describe("cross-PCF on a domain", () => {
    const n = 60_000;
    const x = new Float32Array(n);
    const y = new Float32Array(n);
    const type = new Uint8Array(n);
    const region = new Uint8Array(n);
    let seed = 1;
    const rand = () => (seed = (seed * 16807) % 2147483647) / 2147483647;
    for (let i = 0; i < n; i++) {
        x[i] = rand() * 600 - 100; // some cells outside the domain
        y[i] = rand() * 600 - 100;
        type[i] = i % 3;
        region[i] = i % 2;
    }
    const base = {
        x: { buffer: x.buffer, datatype: "double" },
        y: { buffer: y.buffer, datatype: "double" },
        type: { buffer: type.buffer, datatype: "text", nTypes: 3 },
        region: { buffer: region.buffer, datatype: "text", code: 1 },
    };

    test("cells go through the same frame as the domain, into typed arrays", () => {
        const domain = polygonGeometryDomain(holed, aniso);
        const cells = gatherCells({ ...base, frame: aniso, domain });
        expect(cells.xs).toBeInstanceOf(Float64Array);
        expect(cells.typeId).toBeInstanceOf(Uint16Array);
        // Every gathered cell is the data-frame cell mapped through the frame.
        const inData = polygonGeometryDomain(holed, identityFrame());
        let insideData = 0;
        for (let i = 0; i < n; i++) if (region[i] === 1 && domainContains(inData, x[i], y[i])) insideData++;
        let insideFrame = 0;
        for (let i = 0; i < cells.xs.length; i++) if (domainContains(domain, cells.xs[i], cells.ys[i])) insideFrame++;
        expect(insideFrame).toBe(insideData);
    });

    test("uniform random cells give g ≈ 1 inside a holed domain under an anisotropic frame", () => {
        // Anisotropy shrinks y 2× relative to x: the domain and the cells must both see it, or the
        // density is wrong near the hole and g drifts from 1.
        const res = crossPcfInDomain({
            ...base,
            frame: aniso,
            domain: polygonGeometryDomain(holed, aniso),
            rMax: 10,
            nBins: 4,
            edgeCorrected: true,
        });
        const inside = res.counts.reduce((a, b) => a + b, 0);
        expect(inside / (n / 2)).toBeCloseTo((400 ** 2 - 100 ** 2) / 600 ** 2, 1);
        expect(res.candidates - (res.outside ?? 0)).toBe(inside);
        let worst = 0;
        for (const v of res.g) worst = Math.max(worst, Math.abs(v - 1));
        expect(worst).toBeLessThan(0.15);
    });
});

describe("physical unit override", () => {
    test("scales after the store's transform and says so", () => {
        const f = withPhysicalUnit(frameFromMatrix(new Matrix4().translate([100, 0, 0]), "px", "store"), 0.5, "test");
        expect(f.unit).toBe("µm");
        expect(f.provenance).toContain("test");
        expect(f.toAnalysis[0]).toBeCloseTo(0.5);
        expect(f.toAnalysis[12]).toBeCloseTo(50);
    });
});
