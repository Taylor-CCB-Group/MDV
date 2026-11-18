import { describe, expect, test } from "vitest";
import { createRowSelector } from "@/links/link_utils";

/** A promise we resolve by hand, to control the order in which async filtered-indices requests finish. */
function deferred<T>() {
    let resolve!: (v: T) => void;
    const promise = new Promise<T>((r) => {
        resolve = r;
    });
    return { promise, resolve };
}

/** `filterArray[i] === 0` means row i is visible, as in DataStore. */
function makeSource(opts: { size?: number; visible: number[]; highlighted?: number[] }) {
    const filterArray = new Uint8Array(opts.size ?? 10).fill(1);
    for (const i of opts.visible) filterArray[i] = 0;
    const requests: ReturnType<typeof deferred<Uint32Array>>[] = [];
    const source = {
        filterArray,
        highlighted: opts.highlighted,
        getHighlightedData() {
            return source.highlighted;
        },
        getFilteredIndices() {
            const d = deferred<Uint32Array>();
            requests.push(d);
            return d.promise;
        },
    };
    return { source, requests };
}

function recorder() {
    const applied: number[][] = [];
    return { applied, apply: (rows: ArrayLike<number>) => applied.push(Array.from(rows)) };
}

describe("createRowSelector", () => {
    test("uses the filtered rows when nothing is highlighted", async () => {
        const { source, requests } = makeSource({ visible: [2, 3, 4] });
        const { applied, apply } = recorder();
        const p = createRowSelector(source, apply)();
        requests[0].resolve(new Uint32Array([2, 3, 4]));
        await p;
        expect(applied).toEqual([[2, 3, 4]]);
    });

    test("a highlight takes precedence over the filtered rows, without waiting for the worker", async () => {
        const { source, requests } = makeSource({ visible: [2, 3, 4], highlighted: [3] });
        const { applied, apply } = recorder();
        await createRowSelector(source, apply)();
        expect(applied).toEqual([[3]]);
        expect(requests).toHaveLength(0);
    });

    test("only the highlighted rows that pass the filter are used", async () => {
        const { source } = makeSource({ visible: [2, 3, 4], highlighted: [1, 3, 9] });
        const { applied, apply } = recorder();
        await createRowSelector(source, apply)();
        expect(applied).toEqual([[3]]);
    });

    test("falls back to the filtered rows when the highlight is entirely filtered out", async () => {
        const { source, requests } = makeSource({ visible: [2, 3], highlighted: [7] });
        const { applied, apply } = recorder();
        const p = createRowSelector(source, apply)();
        requests[0].resolve(new Uint32Array([2, 3]));
        await p;
        expect(applied).toEqual([[2, 3]]);
    });

    test("a stale filtered result that finishes late does not overwrite a newer highlight", async () => {
        const { source, requests } = makeSource({ visible: [2, 3, 4] });
        const { applied, apply } = recorder();
        const refresh = createRowSelector(source, apply);
        const first = refresh(); // 'filtered' - waiting on the worker
        source.highlighted = [3];
        const second = refresh(); // 'data_highlighted' arrives before the worker replies
        requests[0].resolve(new Uint32Array([2, 3, 4]));
        await Promise.all([first, second]);
        expect(applied).toEqual([[3]]);
    });

    test("filtered results finishing out of order: the latest request wins", async () => {
        const { source, requests } = makeSource({ visible: [2, 3, 4] });
        const { applied, apply } = recorder();
        const refresh = createRowSelector(source, apply);
        const first = refresh();
        const second = refresh();
        requests[1].resolve(new Uint32Array([4]));
        requests[0].resolve(new Uint32Array([2, 3, 4]));
        await Promise.all([first, second]);
        expect(applied).toEqual([[4]]);
    });

    test("a superseded refresh resolves only after the latest has been applied", async () => {
        const { source, requests } = makeSource({ visible: [2, 3, 4] });
        const { applied, apply } = recorder();
        const refresh = createRowSelector(source, apply);
        let firstDone = false;
        const first = refresh().then(() => {
            firstDone = true;
        });
        const second = refresh();
        requests[0].resolve(new Uint32Array([2, 3, 4]));
        await new Promise((r) => setTimeout(r, 0));
        // `first` was superseded, but must not resolve while the result it deferred to is still outstanding
        expect(firstDone).toBe(false);
        expect(applied).toEqual([]);
        requests[1].resolve(new Uint32Array([4]));
        await Promise.all([first, second]);
        expect(firstDone).toBe(true);
        expect(applied).toEqual([[4]]);
    });

    test("clearing the highlight goes back to the filtered rows", async () => {
        const { source, requests } = makeSource({ visible: [2, 3], highlighted: [3] });
        const { applied, apply } = recorder();
        const refresh = createRowSelector(source, apply);
        await refresh();
        source.highlighted = [];
        const p = refresh();
        requests[0].resolve(new Uint32Array([2, 3]));
        await p;
        expect(applied).toEqual([[3], [2, 3]]);
    });
});
