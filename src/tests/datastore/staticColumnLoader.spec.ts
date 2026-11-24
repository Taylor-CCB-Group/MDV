import pako from "pako";
import { afterEach, describe, expect, test, vi } from "vitest";
import type { DataColumn, DataType } from "@/charts/charts";
import { getLocalCompressedBinaryDataLoader } from "@/dataloaders/DataLoaders";

const size = 3;
const good = new Float32Array([1, 2, 3]);
const goodBytes = pako.deflate(new Uint8Array(good.buffer));
const badBytes = new Uint8Array([1, 2, 3, 4, 5]); // not zlib
const body = new Uint8Array(goodBytes.length + badBytes.length);
body.set(goodBytes, 0);
body.set(badBytes, goodBytes.length);
const index = {
    good: [0, goodBytes.length - 1],
    bad: [goodBytes.length, body.length - 1],
};

function column(field: string): DataColumn<DataType> {
    return { field, name: field, datatype: "double" } as DataColumn<DataType>;
}

function stubStaticServer() {
    vi.stubGlobal(
        "fetch",
        vi.fn(async (url: string, args?: { headers?: { range?: string } }) => {
            if (url.endsWith(".json")) return new Response(JSON.stringify(index));
            const [start, end] = (args?.headers?.range ?? "").replace("bytes=", "").split("-").map(Number);
            return new Response(body.slice(start, end + 1));
        }),
    );
}

describe("static column loader", () => {
    afterEach(() => {
        vi.unstubAllGlobals();
    });

    // ChartManager retries any requested column a loader leaves out
    test("leaves out a column it cannot inflate and returns the rest", async () => {
        stubStaticServer();
        const loader = getLocalCompressedBinaryDataLoader([{ name: "cells", size } as never], "/static");
        vi.spyOn(console, "warn").mockImplementation(() => {});

        const result = await loader([column("good"), column("bad")], "cells", size);

        expect(result.map((c) => c.field)).toEqual(["good"]);
        const data = result[0].data;
        expect(data).toBeInstanceOf(SharedArrayBuffer);
        expect(Array.from(new Float32Array(data as SharedArrayBuffer))).toEqual([1, 2, 3]);
    });
});
