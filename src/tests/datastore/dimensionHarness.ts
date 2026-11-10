import { expect, vi } from "vitest";
import DataStore from "@/datastore/DataStore";
import type Dimension from "@/datastore/Dimension";
import {
    type FilteredIndexRequest,
    type FilteredIndexResponse,
    handleFilteredIndexRequest,
} from "@/datastore/filteredIndexEval";
import { type FilterRequest, type FilterResponse, handleFilterRequest } from "@/datastore/rangeFilterEval";

/**
 * Shared helpers for Dimension filtering tests.
 *
 * Dimension subclasses and DataStore construct Workers eagerly; jsdom has none, so callers must
 * `stubWorker()` before constructing a store (the aggregation workers are not under test here).
 *
 * With `rangeFilter: true`, the range filter worker is replaced by {@link InProcessFilterWorker},
 * which runs the real message handler a task later, like a worker would.
 *
 * With `filteredIndex: true`, the store's filtered index worker is replaced by
 * {@link InProcessFilteredIndexWorker}, which runs the real message handler when the test tells it to.
 */
export function stubWorker({ rangeFilter = false, filteredIndex = false } = {}) {
    const rangeFilterWorkers: InProcessFilterWorker[] = [];
    const filteredIndexWorkers: InProcessFilteredIndexWorker[] = [];
    vi.stubGlobal(
        "Worker",
        vi.fn().mockImplementation(function MockWorker(url: URL | string) {
            if (filteredIndex && String(url).includes("filteredIndexWorker")) {
                const worker = new InProcessFilteredIndexWorker();
                filteredIndexWorkers.push(worker);
                return worker;
            }
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
    return { rangeFilterWorkers, filteredIndexWorkers };
}

/**
 * Stands in for `filteredIndexWorker.ts`. Like a worker, it handles requests one at a time in the order they were
 * posted, reading the shared filter as it is at that moment - but only when the test calls {@link runNext}, so that
 * a test can change the filter between a request being posted and the worker getting to it.
 */
export class InProcessFilteredIndexWorker {
    onmessage: ((event: { data: FilteredIndexResponse }) => void) | null = null;
    onerror: ((event: unknown) => void) | null = null;
    onmessageerror: ((event: unknown) => void) | null = null;
    /** posted and not yet handled */
    queue: FilteredIndexRequest[] = [];
    posted = 0;
    terminated = false;

    postMessage(message: FilteredIndexRequest) {
        this.posted++;
        this.queue.push(structuredClone(message));
    }

    /** Handles the oldest waiting request and delivers the reply. */
    runNext() {
        const request = this.queue.shift();
        if (!request) throw new Error("no request waiting");
        this.onmessage?.({ data: handleFilteredIndexRequest(request) });
    }

    /** The worker dies (e.g. its script failed to load) without answering what is waiting. */
    fail() {
        this.queue.length = 0;
        this.onerror?.(new Error("worker failed"));
    }

    terminate() {
        this.terminated = true;
    }
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
            throw new Error(`row ${i}: parent count ${ds.filterArray[i]} but ${excluding} dimension(s) exclude it`);
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
