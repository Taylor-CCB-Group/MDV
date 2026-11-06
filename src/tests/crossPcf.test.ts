import type { SimpleFeatureCollection } from "@deck.gl-community/editable-layers";
import { Matrix4 } from "@math.gl/core";
import { domainArea, domainContains } from "intraspatial/spatial/domain";
import { describe, expect, test } from "vitest";
import type { Gate } from "@/react/gates/types";
import {
    type AnalysisFrame,
    frameFromMatrix,
    identityFrame,
    regionsScaleFrame,
    withPhysicalUnit,
} from "@/spatial_stats/analysisFrame";
import { matrixAt, pairCurve, pickPair, resolvePairs, togglePair } from "@/spatial_stats/CrossPcfMatrix";
import { crossPcfInDomain, gatherCells } from "@/spatial_stats/crossPcfCompute";
import type { DomainSource } from "@/spatial_stats/domainSource";
import {
    featureRings,
    gateDomainSource,
    polygonGeometryDomain,
    selectionDomainSource,
} from "@/spatial_stats/domainSource";
import {
    defaultDomainSourceRef,
    listDomainSources,
    refKey,
    resolveDomainSourceRef,
} from "@/spatial_stats/domainSourceRef";
import { boundsDomain, imageExtentFromShape } from "@/spatial_stats/regionDomain";
import { regionImageKey, type SpatialDataTransforms, spatialDataRegionFrame } from "@/spatial_stats/spatialDataFrame";
import { tweenValue } from "@/spatial_stats/useTweened";

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

/** Whole-region stand-in: a fixed square per region, labelled so tests can see which one was used. */
function fakeRegionSource(region?: string): DomainSource {
    const bounds = { minX: 0, minY: 0, maxX: 500, maxY: 500 };
    return {
        label: `whole ${region ?? "extent"}`,
        columns: ["x", "y"],
        region,
        toDomain: (f) => boundsDomain(bounds, f),
    };
}
function regionEnv(regions: string[] = ["r1"]) {
    const made = new Map<string, DomainSource>();
    return {
        regions,
        regionSource: (region?: string) => {
            const k = region ?? "";
            if (!made.has(k)) made.set(k, fakeRegionSource(region));
            return made.get(k) ?? fakeRegionSource(region);
        },
    };
}

describe("domain source refs", () => {
    const ds = {};
    const other = {};
    const empty: SimpleFeatureCollection = { type: "FeatureCollection", features: [] };
    const scatter = { id: "s1", title: "Spatial", param: ["x", "y"], region: "r1", selectionFeatureCollection: holed };
    const env = (selection = holed) => ({
        dataStore: ds,
        charts: [
            { dataStore: ds, config: { ...scatter, selectionFeatureCollection: selection } },
            { dataStore: ds, config: { id: "t", param: ["a"] } }, // no selection layer
            { dataStore: other, config: { ...scatter, id: "s2" } }, // other datasource
        ],
        gates: new Map<string, Gate>([
            ["g1", { id: "g1", name: "G", geometry: holed, columns: ["x", "y"], createdAt: 0, labelPosition: [0, 0] }],
        ]),
        ...regionEnv(),
    });

    test("lists selection hosts on the same datastore, its gates and its regions", () => {
        const list = listDomainSources(env());
        expect(list.map((s) => refKey(s.ref))).toEqual(["selection:s1", "gate:g1", "region:r1"]);
        expect(listDomainSources(env(empty))[0].ready).toBe(false);
        expect(listDomainSources({ ...env(), regions: [] }).at(-1)?.label).toBe("Whole extent");
    });

    test("default is a drawn selection, else the whole region", () => {
        expect(defaultDomainSourceRef(env())).toEqual({ kind: "selection", chartId: "s1" });
        expect(defaultDomainSourceRef(env(empty))).toEqual({ kind: "region", region: "r1" });
        expect(defaultDomainSourceRef({ ...env(empty), regions: [] })).toEqual({ kind: "region" });
    });

    test("resolves to the live geometry; nothing drawn or a deleted gate falls back to the whole region", () => {
        expect(resolveDomainSourceRef({ kind: "selection", chartId: "s1" }, env()).region).toBe("r1");
        expect(resolveDomainSourceRef({ kind: "selection", chartId: "s1" }, env(empty)).label).toBe("whole r1");
        expect(resolveDomainSourceRef({ kind: "gate", gateId: "nope" }, env()).label).toBe("whole extent");
        expect(resolveDomainSourceRef({ kind: "gate", gateId: "g1" }, env()).columns).toEqual(["x", "y"]);
        expect(resolveDomainSourceRef({ kind: "region", region: "r1" }, env()).label).toBe("whole r1");
    });
});

test("a resolved source is stable for the same geometry and new for an edit", () => {
    const config = { id: "s", param: ["x", "y"], selectionFeatureCollection: holed };
    const env = { dataStore: {}, charts: [{ dataStore: {}, config }], gates: new Map<string, Gate>(), ...regionEnv() };
    const ref = { kind: "selection", chartId: "s" } as const;
    const a = resolveDomainSourceRef(ref, env);
    expect(resolveDomainSourceRef(ref, env)).toBe(a);
    config.selectionFeatureCollection = { ...holed };
    expect(resolveDomainSourceRef(ref, env)).not.toBe(a);
    const r = resolveDomainSourceRef({ kind: "region", region: "r1" }, env);
    expect(resolveDomainSourceRef({ kind: "region", region: "r1" }, env)).toBe(r);
});

describe("whole-region domain", () => {
    test("image extent from an NGFF shape: x last, y second last", () => {
        expect(imageExtentFromShape([1, 13770, 34155])).toEqual({ minX: 0, minY: 0, maxX: 34155, maxY: 13770 });
        expect(imageExtentFromShape([5])).toBeUndefined();
    });

    test("bounds through an axis-aligned frame stay a rectangle in the analysis frame", () => {
        const d = boundsDomain({ minX: 0, minY: 0, maxX: 400, maxY: 400 }, aniso);
        expect(d.kind).toBe("rect");
        expect(d.bbox).toEqual([10, -5, 210, 95]);
    });

    test("bounds through a rotation become the exact quadrilateral", () => {
        const rot = frameFromMatrix(new Matrix4().rotateZ(Math.PI / 4), "µm", "test");
        const d = boundsDomain({ minX: 0, minY: 0, maxX: 100, maxY: 100 }, rot);
        expect(d.kind).toBe("polygon");
        expect(domainArea(d)).toBeCloseTo(100 * 100, 6);
        expect(domainContains(d, 0, 70)).toBe(true); // inside the rotated square, outside its pre-image
    });

    test("uniform cells give g ≈ 1 with the whole extent as the window", () => {
        const n = 40_000;
        const x = new Float32Array(n);
        const y = new Float32Array(n);
        const type = new Uint8Array(n);
        let seed = 7;
        const rand = () => (seed = (seed * 16807) % 2147483647) / 2147483647;
        for (let i = 0; i < n; i++) {
            x[i] = rand() * 400;
            y[i] = rand() * 400;
            type[i] = i % 2;
        }
        const res = crossPcfInDomain({
            x: { buffer: x.buffer, datatype: "double" },
            y: { buffer: y.buffer, datatype: "double" },
            type: { buffer: type.buffer, datatype: "text", nTypes: 2 },
            frame: aniso,
            domain: boundsDomain({ minX: 0, minY: 0, maxX: 400, maxY: 400 }, aniso),
            rMax: 10,
            nBins: 4,
            edgeCorrected: true,
        });
        expect(res.counts[0] + res.counts[1]).toBe(n);
        let worst = 0;
        for (const v of res.g) worst = Math.max(worst, Math.abs(v - 1));
        expect(worst).toBeLessThan(0.1);
    });
});

describe("matrix view", () => {
    // K = 3 types, B = 2 bins; type 1 has no cells.
    const res = {
        nTypes: 3,
        nBins: 2,
        counts: [5, 0, 4],
        r: [0.5, 1.5],
        g: Float64Array.from({ length: 18 }, (_, i) => i),
        pairs: new Float64Array(18).fill(1),
    };

    test("keeps present types, at one bin or the r-mean", () => {
        const bin = matrixAt(res, { kind: "bin", k: 1 });
        expect(bin.types).toEqual([0, 2]);
        expect(Array.from(bin.value)).toEqual([1, 5, 13, 17]); // g[a*6 + b*2 + 1]
        expect(Array.from(matrixAt(res, { kind: "mean" }).value)).toEqual([0.5, 4.5, 12.5, 16.5]);
        expect(pairCurve(res, 2, 0)).toEqual([12, 13]);
    });

    test("chosen pairs are labels; ones naming a missing category are skipped but keep their index", () => {
        const labels = ["acinar", "beta", "ductal"];
        const pairs = [
            { a: "ductal", b: "acinar" },
            { a: "gone", b: "acinar" },
            { a: "beta", b: "beta" },
        ];
        expect(resolvePairs(labels, pairs)).toEqual([
            { a: 2, b: 0, index: 0 },
            { a: 1, b: 1, index: 2 },
        ]);
    });

    test("toggling adds A → B, and removes it again; B → A is a different pair", () => {
        const one = togglePair([], "ductal", "acinar");
        expect(one).toEqual([{ a: "ductal", b: "acinar" }]);
        expect(togglePair(one, "acinar", "ductal")).toHaveLength(2);
        expect(togglePair(one, "ductal", "acinar")).toEqual([]);
    });
});

test("a plain pick replaces the pairs; an additive pick toggles one", () => {
    const two = [
        { a: "ductal", b: "acinar" },
        { a: "beta", b: "beta" },
    ];
    expect(pickPair(two, "beta", "ductal", false)).toEqual([{ a: "beta", b: "ductal" }]);
    expect(pickPair(two, "ductal", "acinar", false)).toEqual([{ a: "ductal", b: "acinar" }]);
    expect(pickPair(two, "beta", "beta", true)).toEqual([{ a: "ductal", b: "acinar" }]);
    expect(pickPair(two, "beta", "ductal", true)).toHaveLength(3);
});

test("scale tween eases out and lands on the target", () => {
    expect(tweenValue(2, 10, 0)).toBe(2);
    expect(tweenValue(2, 10, 1)).toBe(10);
    expect(tweenValue(2, 10, 2)).toBe(10);
    const mid = tweenValue(2, 10, 0.5);
    expect(mid).toBeGreaterThan(6); // ease-out: past halfway at half time
    expect(mid).toBeLessThan(10);
});
