import { Matrix4 } from "@math.gl/core";
import * as glMatrix from "gl-matrix";
import { describe, expect, test } from "vitest";
import { evaluateFilterSpec } from "@/datastore/rangeFilterEval";
import {
    fitOrbitView,
    isSameOrbitView,
    type OrbitSelectionView,
    orbitViewFromLegacyCamera,
    orbitViewport,
    projectRowsToScreen,
    screenProjectionMatrix,
} from "@/react/orbitView";

/**
 * Where the classic 3D scatter drew a data point, in pixels from the top left: the projection in
 * `WGL2DI`'s 3D shader, with the camera `Camera.updateCamera` sets up.
 */
function legacyPixel(
    point: [number, number, number],
    { distance, theta, phi }: { distance: number; theta: number; phi: number },
    center: [number, number, number],
    scales: [number, number, number],
    width: number,
    height: number,
) {
    const c = [-center[0] * 2, center[1] * 2, -center[2] * 2];
    const vf = distance * Math.sin(theta) * Math.cos(phi);
    const vr = distance * Math.cos(theta) * Math.cos(phi);
    const vu = distance * Math.sin(phi);
    // front is +z, right +x, up +y
    const eye = [c[0] + vr, c[1] + vu, c[2] + vf];
    const view = glMatrix.mat4.lookAt(
        glMatrix.mat4.create(),
        glMatrix.vec3.fromValues(eye[0], eye[1], eye[2]),
        glMatrix.vec3.fromValues(c[0], c[1], c[2]),
        glMatrix.vec3.fromValues(0, 1, 0),
    );
    const projection = glMatrix.mat4.perspective(glMatrix.mat4.create(), Math.PI / 4, width / height, 1, 100000);
    const p = [2 * -point[0] * scales[0], 2 * point[1] * scales[1], 2 * -point[2] * scales[2], 1] as glMatrix.vec4;
    const clip = glMatrix.vec4.transformMat4(glMatrix.vec4.create(), p, view);
    glMatrix.vec4.transformMat4(clip, clip, projection);
    const [x, y, , w] = clip;
    return [((x / w + 1) / 2) * width, ((1 - y / w) / 2) * height];
}

const cases = [
    { name: "the default view", camera: { distance: 64, theta: -1.038, phi: 0.261 }, center: [0.1, 0.2, 0.9] },
    { name: "looking along x", camera: { distance: 49.6, theta: 0, phi: 0 }, center: [0, 0, 0] },
    { name: "turned and from below", camera: { distance: 8539, theta: 27.3, phi: -0.6 }, center: [3236, 3159, 767] },
] as const;

describe("a classic 3D scatter's camera", () => {
    test.each(cases)("draws each point where the classic chart did: $name", ({ camera, center }) => {
        const [width, height] = [420, 300];
        const scales: [number, number, number] = [0.5, 1, 3];
        const view = orbitViewFromLegacyCamera(camera, [...center], height);
        expect(view).toBeDefined();
        const selectionView: OrbitSelectionView = { ...view!, width, height, orbitAxis: "Y" };
        const modelMatrix = new Matrix4().scale(scales);
        const r = camera.distance / 8;
        const points: [number, number, number][] = [
            [center[0] / scales[0], center[1] / scales[1], center[2] / scales[2]],
            [center[0] / scales[0] + r, center[1] / scales[1] - r / 2, center[2] / scales[2] + r / 3],
            [center[0] / scales[0] - r, center[1] / scales[1] + r, center[2] / scales[2] - r / 4],
        ];
        const [px, py] = projectRowsToScreen(
            [points.map((p) => p[0]), points.map((p) => p[1]), points.map((p) => p[2])],
            modelMatrix,
            selectionView,
        );
        points.forEach((point, i) => {
            const [lx, ly] = legacyPixel(point, camera, [...center], scales, width, height);
            expect(px[i]).toBeCloseTo(lx, 1);
            expect(py[i]).toBeCloseTo(ly, 1);
        });
    });

    test("is left for the chart to fit when it is incomplete", () => {
        expect(orbitViewFromLegacyCamera({ distance: 0, theta: 0, phi: 0 }, [0, 0, 0], 300)).toBeUndefined();
        expect(orbitViewFromLegacyCamera({ distance: 10, theta: Number.NaN, phi: 0 }, [0, 0, 0], 300)).toBeUndefined();
        expect(orbitViewFromLegacyCamera({ distance: 10, theta: 0, phi: 0 }, [0, 0, 0], 0)).toBeUndefined();
    });
});

test("projecting rows matches deck's own projection", () => {
    const view: OrbitSelectionView = {
        target: [5, -2, 1],
        zoom: 3.2,
        rotationOrbit: 40,
        rotationX: 20,
        width: 500,
        height: 320,
        orbitAxis: "Z",
    };
    const modelMatrix = new Matrix4().scale([2, 1, 0.5]);
    const points = [
        [5, -2, 1],
        [6, -1, 3],
        [3.5, -4, 0],
    ];
    const [px, py] = projectRowsToScreen(
        [points.map((p) => p[0]), points.map((p) => p[1]), points.map((p) => p[2])],
        modelMatrix,
        view,
    );
    const viewport = orbitViewport(view);
    points.forEach((p, i) => {
        const [x, y] = viewport.project(Array.from(modelMatrix.transformAsPoint(p)));
        expect(px[i]).toBeCloseTo(x, 3);
        expect(py[i]).toBeCloseTo(y, 3);
    });
});

test("a fitted view keeps the whole box in the chart however it turns", () => {
    const view = fitOrbitView({ min: [0, 0, 0], max: [10, 20, 5] }, 400, 300);
    expect(view.target).toEqual([5, 10, 2.5]);
    // the box's diagonal fits the shorter side
    expect(Math.hypot(10, 20, 5) * 2 ** view.zoom).toBeLessThan(300);
});

test("a fitted view frames a very small cloud, and a single point", () => {
    const tiny = fitOrbitView({ min: [0, 0, 0], max: [1e-12, 1e-12, 1e-12] }, 300, 300);
    const diagonalPixels = Math.hypot(1e-12, 1e-12, 1e-12) * 2 ** tiny.zoom;
    expect(diagonalPixels).toBeGreaterThan(150);
    expect(diagonalPixels).toBeLessThan(300);
    const point = fitOrbitView({ min: [5, 5, 5], max: [5, 5, 5] }, 300, 300);
    expect(Number.isFinite(point.zoom)).toBe(true);
    expect(point.target).toEqual([5, 5, 5]);
});

describe("the view a selection was drawn in", () => {
    const view: OrbitSelectionView = {
        target: [1, 2, 3],
        zoom: 2,
        rotationOrbit: 10,
        rotationX: 5,
        width: 300,
        height: 200,
        orbitAxis: "Y",
    };

    test("is the same after saving and loading", () => {
        expect(isSameOrbitView(view, JSON.parse(JSON.stringify(view)))).toBe(true);
        expect(isSameOrbitView(view, { ...view, width: 300.6 })).toBe(true);
    });

    test("is not the same once the camera has moved", () => {
        expect(isSameOrbitView(view, { ...view, rotationOrbit: 11 })).toBe(false);
        expect(isSameOrbitView(view, { ...view, zoom: 2.5 })).toBe(false);
        expect(isSameOrbitView(view, { ...view, target: [1, 2.1, 3] })).toBe(false);
        expect(isSameOrbitView(view, { ...view, orbitAxis: "Z" })).toBe(false);
    });
});

test("the filter worker selects what is inside an outline on screen", () => {
    const view: OrbitSelectionView = {
        target: [0, 0, 0],
        zoom: 3,
        rotationOrbit: 30,
        rotationX: 15,
        width: 400,
        height: 300,
        orbitAxis: "Y",
    };
    const modelMatrix = new Matrix4().scale([1, 2, 0.5]);
    // a deterministic spread of points, with a missing value and some behind the camera
    const n = 2000;
    const columns = [0, 1, 2].map((axis) =>
        Float32Array.from({ length: n }, (_, i) => Math.sin(i * (axis + 1) * 12.9898) * 20),
    );
    columns[0][7] = Number.NaN;
    columns[2][11] = 1e6;
    const points: [number, number][] = [
        [120, 80],
        [300, 90],
        [260, 240],
        [140, 200],
        [120, 80],
    ];
    const onWorker = new Uint8Array(n);
    evaluateFilterSpec(
        {
            kind: "screenPoly",
            data1: columns[0],
            data2: columns[1],
            data3: columns[2],
            matrix: Array.from(screenProjectionMatrix(view, modelMatrix)),
            points,
        },
        n,
        onWorker,
    );
    const [px, py] = projectRowsToScreen([columns[0], columns[1], columns[2]], modelMatrix, view);
    const expected = new Uint8Array(n);
    evaluateFilterSpec({ kind: "poly", data1: px, data2: py, points }, n, expected);
    expect(onWorker).toEqual(expected);
    expect(onWorker[7]).toBe(0);
    const selected = onWorker.reduce((a, b) => a + b, 0);
    expect(selected).toBeGreaterThan(50);
    expect(selected).toBeLessThan(n);
});
