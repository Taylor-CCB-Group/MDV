import { beforeEach, describe, expect, test, vi } from "vitest";

// DataStore imports helpers via @/lib/utils
vi.mock("@/lib/utils", async (importOriginal) => {
	const actual: any = await importOriginal();
	return {
		...actual,
		isDatatypeNumeric: (t: string) => !!t.match(/double|float|int/),
	};
});

beforeEach(() => {
	// DataStore uses a WebWorker in its constructor; Node test env doesn't provide it.
	// A minimal stub is enough for these minMax/quantiles tests.
	vi.stubGlobal(
		"Worker",
		class WorkerStub {
			onmessage = null;
			onerror = null;
			postMessage(_msg: any) {}
			terminate() {}
			addEventListener(_type: any, _listener: any) {}
			removeEventListener(_type: any, _listener: any) {}
		},
	);
});

function f32Buffer(values: number[]): ArrayBuffer {
	const buf = new ArrayBuffer(values.length * 4);
	new Float32Array(buf).set(values);
	return buf;
}

describe("DataStore.setColumnData (numeric stats)", () => {
	test("recomputes minMax/quantiles by default when overwriting numeric data", async () => {
		// Import after stubbing Worker
		const { default: DataStore } = await import("../datastore/DataStore.js");
		const ds = new DataStore(3, {
			columns: [{ name: "x", field: "x", datatype: "double" }],
		});

		ds.setColumnData("x", f32Buffer([10, 11, 12]));
		expect(ds.getMinMaxForColumn("x")).toEqual([10, 12]);

		// Overwrite with new data on the same column (DGE-style behavior)
		ds.setColumnData("x", f32Buffer([-5, 0, 5]));
		expect(ds.getMinMaxForColumn("x")).toEqual([-5, 5]);
		// Quantiles should also have been recomputed (existence check)
		expect(ds.columnIndex.x.quantiles).toBeTruthy();
	});

	test("can preserve existing stats when recomputeStats=false", async () => {
		// Import after stubbing Worker
		const { default: DataStore } = await import("../datastore/DataStore.js");
		const ds = new DataStore(3, {
			columns: [{ name: "x", field: "x", datatype: "double" }],
		});

		ds.setColumnData("x", f32Buffer([1, 2, 3]));
		expect(ds.getMinMaxForColumn("x")).toEqual([1, 3]);

		// Force a stale minMax to simulate a caller intentionally setting it
		ds.columnIndex.x.minMax = [100, 200];
		ds.setColumnData("x", f32Buffer([7, 8, 9]), { recomputeStats: false });

		// With recomputeStats disabled, we preserve the existing minMax cache.
		expect(ds.getMinMaxForColumn("x")).toEqual([100, 200]);
	});
});

