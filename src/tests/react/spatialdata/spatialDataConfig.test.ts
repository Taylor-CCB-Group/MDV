import { describe, expect, test } from "vitest";

import { toSerializableSpatialDataViewState } from "@/react/spatialdata/spatialdata_config";

describe("toSerializableSpatialDataViewState", () => {
    test("extracts target + zoom from a view state", () => {
        expect(
            toSerializableSpatialDataViewState({ target: [1, 2, 0], zoom: 3 }),
        ).toEqual({ target: [1, 2, 0], zoom: 3 });
    });

    test("returns null when there is no target", () => {
        expect(toSerializableSpatialDataViewState(null)).toBeNull();
        expect(toSerializableSpatialDataViewState(undefined)).toBeNull();
    });
});
