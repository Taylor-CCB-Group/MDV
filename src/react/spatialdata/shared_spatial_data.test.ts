import { describe, expect, it, vi } from "vitest";

vi.mock("@/dataloaders/DataLoaderUtil", () => ({
    getProjectURL: (path: string) => `http://project.test/${path}/`,
}));

import { spatialDataUrlForRegion } from "./shared_spatial_data";

describe("spatialDataUrlForRegion", () => {
    it("serves the store from the project's spatial folder by default", () => {
        expect(spatialDataUrlForRegion({ spatial: { file: "kidney.zarr" } })).toBe(
            "http://project.test/spatial/kidney.zarr/",
        );
    });

    it("uses spatial.url when set", () => {
        const region = { spatial: { file: "kidney.zarr", url: "http://127.0.0.1:8091/" } };
        expect(spatialDataUrlForRegion(region)).toBe("http://127.0.0.1:8091/");
    });

    it("ignores an empty spatial.url", () => {
        expect(spatialDataUrlForRegion({ spatial: { file: "kidney.zarr", url: "" } })).toBe(
            "http://project.test/spatial/kidney.zarr/",
        );
    });

    it("still needs spatial.file, which names the object", () => {
        expect(spatialDataUrlForRegion({ spatial: { url: "http://127.0.0.1:8091/" } })).toBeNull();
        expect(spatialDataUrlForRegion({})).toBeNull();
    });
});
