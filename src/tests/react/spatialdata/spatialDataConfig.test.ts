import {
    RENDER_STACK_SCHEMA_VERSION,
    type RenderStack,
} from "@spatialdata/layers";
import { describe, expect, test } from "vitest";

import {
    removeSpatialDataRootViv,
    toSerializableSpatialDataViewState,
} from "@/react/spatialdata/spatialdata_config";

function imageStackWithChannels(): RenderStack {
    return {
        schemaVersion: RENDER_STACK_SCHEMA_VERSION,
        entries: [
            {
                kind: "spatial",
                id: "spatialdata-image-imageA",
                visible: true,
                source: { elementType: "image", elementKey: "imageA" },
                props: {
                    opacity: 1,
                    channels: {
                        channelsVisible: [true, false],
                        colors: [[255, 0, 0], [0, 255, 0]],
                    },
                },
            },
        ],
    };
}

describe("SpatialData chart config ownership", () => {
    test("uses root viewState and layer channels without root viv state", () => {
        const renderStack = imageStackWithChannels();
        const config = {
            viewState: toSerializableSpatialDataViewState({
                target: [1, 2, 0],
                zoom: 3,
            }),
            renderStack,
            viv: {
                channelsStore: { channelsVisible: [false] },
                viewerStore: { viewState: { target: [999, 999, 0], zoom: -1 } },
            },
        };

        const serialized = removeSpatialDataRootViv(config);

        expect("viv" in serialized).toBe(false);
        expect(serialized.viewState).toEqual({ target: [1, 2, 0], zoom: 3 });
        expect(serialized.renderStack.entries[0]?.props.channels).toEqual({
            channelsVisible: [true, false],
            colors: [[255, 0, 0], [0, 255, 0]],
        });
    });
});
