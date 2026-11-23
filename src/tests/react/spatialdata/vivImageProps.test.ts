import { RENDER_STACK_SCHEMA_VERSION, type RenderStack } from "@spatialdata/layers";
import { describe, expect, test } from "vitest";
import { vivImageLayerProps } from "@/react/spatialdata/viv_image_props";

function stackWithImage(props: Record<string, unknown>): RenderStack {
    return {
        schemaVersion: RENDER_STACK_SCHEMA_VERSION,
        entries: [
            {
                kind: "spatial",
                id: "image-1",
                visible: true,
                source: { elementType: "image", elementKey: "img" },
                props: { opacity: 1, ...props },
            },
            { kind: "host", id: "deck:scatter", visible: true, source: { hostLayerId: "scatter" }, props: {} },
        ],
    };
}

describe("vivImageLayerProps", () => {
    // an image pick stole clicks meant for the points under it (on WebGPU, where Viv's tiles pick)
    test("an image is never pickable", () => {
        expect(vivImageLayerProps(undefined, "image-1")).toEqual({ pickable: false });
        expect(vivImageLayerProps(stackWithImage({}), "image-1")).toEqual({ pickable: false });
        expect(vivImageLayerProps(stackWithImage({ vivLayerProps: { pickable: true } }), "image-1")).toEqual({
            pickable: false,
        });
    });

    test("keeps the props saved on the entry", () => {
        const vivLayerProps = { brightness: [0.4], contrast: [0.6] };
        const stack = stackWithImage({ vivLayerProps });
        expect(vivImageLayerProps(stack, "image-1")).toEqual({ ...vivLayerProps, pickable: false });
        // a copy: the entry's props are patched in place, and must not gain `pickable`
        expect(vivLayerProps).not.toHaveProperty("pickable");
    });

    test("ignores entries that are not spatial", () => {
        expect(vivImageLayerProps(stackWithImage({}), "deck:scatter")).toEqual({ pickable: false });
    });
});
