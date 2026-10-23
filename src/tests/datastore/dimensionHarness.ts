import { expect, vi } from "vitest";
import DataStore from "@/datastore/DataStore";
import type Dimension from "@/datastore/Dimension";
import { type FilterRequest, type FilterResponse, handleFilterRequest } from "@/datastore/rangeFilterEval";

/**
 * Shared helpers for Dimension filtering tests.
 *
 * Dimension subclasses and DataStore construct Workers eagerly; jsdom has none, so callers must
 * `stubWorker()` before constructing a store (the aggregation workers are not under test here).
 *
 * With `rangeFilter: true`, the range filter worker is replaced by {@link InProcessFilterWorker},
 * which runs the real message handler a task later, like a worker would.
 */
export function stubWorker({ rangeFilter = false } = {}) {
    const rangeFilterWorkers: InProcessFilterWorker[] = [];
    vi.stubGlobal(
        "Worker",
        vi.fn().mockImplementation(function MockWorker(url: URL | string) {
            if (rangeFilter && String(url).includes("rangeFilterWorker")) {
                const worker = new InProcessFilterWorker();
                rangeFilterWorkers.push(worker);
                return worker;
            }
            return {
                addEventListener: vi.fn(),
                removeEventListener: vi.fn(),
                postMessage: vi.fn(),
                terminate: vi.fn(),
            };
        }),
    );
    return { rangeFilterWorkers };
}

/**
 * Stands in for `rangeFilterWorker.ts`: clones the request as `postMessage` would (so shared
 * buffers stay shared and unclonable specs throw), and responds a task later unless terminated.
 */
export class InProcessFilterWorker {
    onmessage: ((event: { data: FilterResponse }) => void) | null = null;
    onerror: ((event: unknown) => void) | null = null;
    requests: FilterRequest[] = [];
    terminated = false;

    postMessage(message: FilterRequest) {
        const request = structuredClone(message);
        this.requests.push(request);
        setTimeout(() => {
            if (!this.terminated) this.onmessage?.({ data: handleFilterRequest(request) });
        }, 0);
    }

    terminate() {
        this.terminated = true;
    }
}

/** Resolves after every dimension's pending async filter (and anything it queued) has landed. */
export async function settleFilters(ds: DataStore) {
    do {
        await new Promise((resolve) => setTimeout(resolve, 0));
    } while ((ds.dimensions as Dimension[]).some((d) => d.pendingFilter));
}

export type ColumnSpec =
    | { field: string; datatype: "double" | "integer" | "int32"; data: number[] }
    | { field: string; datatype: "text" | "text16"; data: string[] };

/** Builds a real DataStore with the given columns loaded. */
export function createStore(size: number, columns: ColumnSpec[]) {
    const ds = new DataStore(size, {
        name: "test",
        columns: columns.map((c) => ({
            field: c.field,
            name: c.field,
            datatype: c.datatype,
        })),
    } as never);
    for (const c of columns) {
        ds.setColumnData(c.field, c.data);
    }
    return ds;
}

/** Rows that pass every dimension's filter (global view). */
export function passingRows(ds: DataStore) {
    const rows: number[] = [];
    for (let i = 0; i < ds.size; i++) {
        if (ds.filterArray[i] === 0) rows.push(i);
    }
    return rows;
}

/**
 * Checks the DataStore/Dimension bookkeeping against a reference model:
 * - every local state is one of 0/1/2/3, and 2/3 only appear with a background filter;
 * - `parent.filterArray[i]` is the number of dimensions excluding row i (local state 1 or 3);
 * - `parent.filterSize` is the number of rows nobody excludes.
 */
export function assertConsistent(ds: DataStore) {
    const dims = ds.dimensions as Dimension[];
    let passing = 0;
    for (let i = 0; i < ds.size; i++) {
        let excluding = 0;
        for (const d of dims) {
            const s = d.filterArray[i];
            if (s > 3) throw new Error(`row ${i}: invalid local state ${s}`);
            if (s >= 2 && !d.bgfData) {
                throw new Error(`row ${i}: background state ${s} without background filter`);
            }
            if (s === 1 || s === 3) excluding++;
        }
        if (ds.filterArray[i] !== excluding) {
            throw new Error(
                `row ${i}: parent count ${ds.filterArray[i]} but ${excluding} dimension(s) exclude it`,
            );
        }
        if (excluding === 0) passing++;
    }
    expect(ds.filterSize).toBe(passing);
}

/** Deterministic PRNG (mulberry32) so randomized tests are reproducible. */
export function rng(seed: number) {
    let a = seed >>> 0;
    return () => {
        a = (a + 0x6d2b79f5) >>> 0;
        let t = a;
        t = Math.imul(t ^ (t >>> 15), t | 1);
        t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
        return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
}
