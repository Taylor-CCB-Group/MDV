import { describe, expect, test, vi, beforeEach, afterEach } from "vitest";
import {
	computeEffectSize,
	log2FoldChange,
	computeGeneStats,
} from "@/datastore/dgeStats";
import { runDGEOnDataStore } from "@/datastore/dgeIntegration";
import * as linkUtils from "@/links/link_utils";
import { DGERunner } from "@/datastore/DGEDimension";

// ── Data type detection logic ───────────────────────────────────────────────
// Mirrors the multi-gene probe loop in dgeIntegration.ts: runDGEOnDataStore.
// The real code scans values across multiple genes, tracking globalMaxVal and
// globalHasNegative, and only makes a decision after examining all probes (or
// early-exits on a definitive signal).

function resolveDataTypeMultiGene(
	geneArrays: (Float32Array | null)[],
): "log1p" | "linear" | "zscored" {
	let globalMaxVal = -Infinity;
	let globalHasNegative = false;
	let anyFinite = false;

	for (const arr of geneArrays) {
		if (!arr) continue;
		for (let j = 0; j < arr.length; j++) {
			const v = arr[j];
			if (Number.isNaN(v)) continue;
			anyFinite = true;
			if (v < 0) { globalHasNegative = true; break; }
			if (v > globalMaxVal) globalMaxVal = v;
		}
		if (globalHasNegative || globalMaxVal > 20) break;
	}

	if (!anyFinite) return "log1p";
	if (globalHasNegative) return "zscored";
	return globalMaxVal > 20 ? "linear" : "log1p";
}

describe("multi-gene data type detection", () => {
	test("detects linear when LATER genes have high values (raw count data)", () => {
		const gene0 = new Float32Array([NaN, NaN, NaN, NaN]);
		const gene1 = new Float32Array([NaN, 1, NaN, 2]);
		const gene2 = new Float32Array([NaN, NaN, 3, NaN]);
		const gene3 = new Float32Array([NaN, 50, NaN, 116]);
		expect(resolveDataTypeMultiGene([gene0, gene1, gene2, gene3])).toBe("linear");
	});

	test("first gene with low values does NOT prematurely decide log1p", () => {
		const gene0 = new Float32Array([NaN, 2.0, NaN, NaN]);
		const gene1 = new Float32Array([NaN, NaN, NaN, NaN]);
		const gene2 = new Float32Array([NaN, 100, NaN, NaN]);
		expect(resolveDataTypeMultiGene([gene0, gene1, gene2])).toBe("linear");
	});

	test("detects log1p when ALL genes have small values", () => {
		const gene0 = new Float32Array([NaN, 0.5, NaN, 1.2]);
		const gene1 = new Float32Array([NaN, NaN, 3.1, NaN]);
		expect(resolveDataTypeMultiGene([gene0, gene1])).toBe("log1p");
	});

	test("detects zscored immediately on negative value", () => {
		const gene0 = new Float32Array([NaN, -0.5, NaN, 2.3]);
		const gene1 = new Float32Array([NaN, 50, NaN, NaN]);
		expect(resolveDataTypeMultiGene([gene0, gene1])).toBe("zscored");
	});

	test("all-NaN across all genes defaults to log1p", () => {
		const gene0 = new Float32Array([NaN, NaN, NaN]);
		const gene1 = new Float32Array([NaN, NaN]);
		expect(resolveDataTypeMultiGene([gene0, gene1])).toBe("log1p");
	});

	test("null buffers are skipped gracefully", () => {
		const gene0 = null;
		const gene1 = new Float32Array([NaN, 25, NaN]);
		expect(resolveDataTypeMultiGene([gene0, gene1])).toBe("linear");
	});
});

describe("computeEffectSize", () => {
	test("uses log2FC formula for log1p data", () => {
		const result = computeEffectSize(2.0, 1.0, "log1p");
		const expected = log2FoldChange(2.0, 1.0, true);
		expect(result).toBe(expected);
	});

	test("uses log2FC formula for linear data (without expm1)", () => {
		const result = computeEffectSize(2.0, 1.0, "linear");
		const expected = log2FoldChange(2.0, 1.0, false);
		expect(result).toBe(expected);
	});

	test("uses mean difference for z-scored data", () => {
		const result = computeEffectSize(2.0, 1.0, "zscored");
		expect(result).toBe(1.0);
	});

	test("mean difference works with negative values (z-scored)", () => {
		const result = computeEffectSize(-0.5, 0.3, "zscored");
		expect(result).toBeCloseTo(-0.8, 10);
	});

	test("mean difference is zero when means are equal", () => {
		expect(computeEffectSize(3.5, 3.5, "zscored")).toBe(0);
	});

	test("log2FC is zero when means are equal", () => {
		expect(computeEffectSize(1.5, 1.5, "log1p")).toBeCloseTo(0, 6);
	});

	test("log2FC is positive when target > reference (log1p data)", () => {
		expect(computeEffectSize(3.0, 1.0, "log1p")).toBeGreaterThan(0);
	});

	test("mean diff is positive when target > reference (z-scored data)", () => {
		expect(computeEffectSize(3.0, 1.0, "zscored")).toBeGreaterThan(0);
	});
});

describe("computeGeneStats with dataType parameter", () => {
	function makeArrays(values: number[], groups: number[]) {
		return {
			v: new Float32Array(values),
			g: new Uint8Array(groups),
			f: new Uint8Array(values.length),
		};
	}

	test("uses log2FC for log1p data type", () => {
		const values = [2, 3, 4, 5, 6, 10, 11, 12, 13, 14];
		const groups = [0, 0, 0, 0, 0, 1, 1, 1, 1, 1];
		const { v, g, f } = makeArrays(values, groups);

		const result = computeGeneStats("Gene1", v, g, f, 0, -1, "log1p");
		const expectedEffectSize = log2FoldChange(4, 12, true);
		expect(result.effectSize).toBeCloseTo(expectedEffectSize, 3);
	});

	test("uses log2FC for linear data type (without expm1)", () => {
		const values = [2, 3, 4, 5, 6, 10, 11, 12, 13, 14];
		const groups = [0, 0, 0, 0, 0, 1, 1, 1, 1, 1];
		const { v, g, f } = makeArrays(values, groups);

		const result = computeGeneStats("Gene1", v, g, f, 0, -1, "linear");
		const expectedEffectSize = log2FoldChange(4, 12, false);
		expect(result.effectSize).toBeCloseTo(expectedEffectSize, 3);
	});

	test("uses mean difference for zscored data type", () => {
		const values = [-1, -0.5, 0, 0.5, 1, 2, 2.5, 3, 3.5, 4];
		const groups = [0, 0, 0, 0, 0, 1, 1, 1, 1, 1];
		const { v, g, f } = makeArrays(values, groups);

		const result = computeGeneStats("Gene1", v, g, f, 0, -1, "zscored");
		expect(result.effectSize).toBeCloseTo(0 - 3, 3);
	});

	test("effect size sign is consistent across data types", () => {
		const values = [1, 1, 1, 1, 1, 5, 5, 5, 5, 5];
		const groups = [0, 0, 0, 0, 0, 1, 1, 1, 1, 1];
		const { v, g, f } = makeArrays(values, groups);

		const resultLog1p = computeGeneStats("Gene1", v, g, f, 0, -1, "log1p");
		const resultLinear = computeGeneStats("Gene1", v, g, f, 0, -1, "linear");
		const resultZscored = computeGeneStats("Gene1", v, g, f, 0, -1, "zscored");

		expect(resultLog1p.effectSize).toBeLessThan(0);
		expect(resultLinear.effectSize).toBeLessThan(0);
		expect(resultZscored.effectSize).toBeLessThan(0);
	});
});

describe("gene ordering alignment", () => {
	test("maps DGE results to correct gene row indices", () => {
		const valueToRowIndex = new Map<string, number>([
			["GeneA", 0],
			["GeneB", 1],
			["GeneC", 2],
			["GeneD", 3],
			["GeneE", 4],
		]);
		const n = 5;

		const dgeResults = [
			{ gene: "GeneD", effectSize: 2.5, pval: 0.001 },
			{ gene: "GeneA", effectSize: -1.0, pval: 0.01 },
			{ gene: "GeneC", effectSize: 0.5, pval: 0.05 },
		];

		const effectSizeColumn = new Float32Array(n).fill(NaN);
		const pvalColumn = new Float32Array(n).fill(NaN);

		for (const result of dgeResults) {
			const idx = valueToRowIndex.get(result.gene);
			if (idx !== undefined) {
				effectSizeColumn[idx] = result.effectSize;
				pvalColumn[idx] = result.pval;
			}
		}

		expect(effectSizeColumn[0]).toBeCloseTo(-1.0, 5);
		expect(pvalColumn[0]).toBeCloseTo(0.01, 5);
		expect(effectSizeColumn[1]).toBeNaN();
		expect(pvalColumn[1]).toBeNaN();
		expect(effectSizeColumn[2]).toBeCloseTo(0.5, 5);
		expect(effectSizeColumn[3]).toBeCloseTo(2.5, 5);
		expect(effectSizeColumn[4]).toBeNaN();
	});

	test("handles duplicate gene names by using last occurrence", () => {
		const valueToRowIndex = new Map<string, number>([
			["Gene1", 0],
			["Gene2", 1],
		]);

		const dgeResults = [
			{ gene: "Gene1", effectSize: 1.0 },
			{ gene: "Gene1", effectSize: 2.0 },
		];

		const column = new Float32Array(2).fill(NaN);
		for (const result of dgeResults) {
			const idx = valueToRowIndex.get(result.gene);
			if (idx !== undefined) {
				column[idx] = result.effectSize;
			}
		}

		expect(column[0]).toBeCloseTo(2.0, 5);
	});

	test("handles genes not in the DataStore gracefully", () => {
		const valueToRowIndex = new Map<string, number>([
			["Gene1", 0],
		]);

		const column = new Float32Array(1).fill(NaN);
		const idx = valueToRowIndex.get("UnknownGene");
		expect(idx).toBeUndefined();
		expect(column[0]).toBeNaN();
	});
});

describe("runDGEOnDataStore link selection", () => {
	beforeEach(() => {
		vi.restoreAllMocks();
	});

	afterEach(() => {
		vi.restoreAllMocks();
	});

	function makeGenesDs(size: number) {
		const columnIndex: Record<string, any> = {};
		return {
			size,
			columnIndex,
			addColumn: vi.fn((meta: any, buffer: SharedArrayBuffer) => {
				columnIndex[meta.field] = { ...meta, buffer, data: new Float32Array(buffer) };
			}),
			setColumnData: vi.fn((field: string, buffer: SharedArrayBuffer) => {
				columnIndex[field] = { ...(columnIndex[field] ?? { field }), buffer, data: new Float32Array(buffer) };
			}),
			dataChanged: vi.fn(),
		};
	}

	test("uses selected genes datasource link when multiple rows_as_columns links exist", async () => {
		const nCells = 3;
		const groupBuffer = new SharedArrayBuffer(nCells);
		new Uint8Array(groupBuffer).set([0, 1, 0]);
		const filterBuffer = new SharedArrayBuffer(nCells);

		const cellsDs: any = {
			size: nCells,
			name: "cells",
			columnIndex: {
				group: {
					datatype: "text",
					values: ["target", "ref"],
					buffer: groupBuffer,
				},
			},
			filterBuffer,
		};

		const genesRna = makeGenesDs(3);
		const genesProtein = makeGenesDs(3);
		const linksByField: Record<string, any> = {
			rna: {
				linkedDs: { name: "rna", dataStore: { columnIndex: { gene_name: { data: new Uint8Array([0, 1, 2]) } } } },
				link: {
					name_column: "gene_name",
					subgroups: { sg_rna: { type: "sparse", name: "RNA", label: "RNA" } },
					valueToRowIndex: new Map<string, number>([["GeneA", 0], ["GeneB", 1]]),
					initPromise: Promise.resolve(),
				},
			},
			protein: {
				linkedDs: { name: "protein", dataStore: { columnIndex: { feature_name: { data: new Uint8Array([0, 1, 2]) } } } },
				link: {
					name_column: "feature_name",
					subgroups: { sg_protein: { type: "sparse", name: "Protein", label: "Protein" } },
					valueToRowIndex: new Map<string, number>([["GeneA", 2], ["GeneB", 0]]),
					initPromise: Promise.resolve(),
				},
			},
		};

		vi.spyOn(linkUtils, "getRowsAsColumnsLinks").mockImplementation(() => [
			linksByField.rna,
			linksByField.protein,
		]);
		vi.spyOn(linkUtils, "getFieldName").mockImplementation((sg, value, idx) => `${sg}|${value}|${idx}`);
		vi.spyOn(DGERunner.prototype, "run").mockResolvedValue({
			results: [
				{
					gene: "GeneB",
					pval: 0.01,
					pvalAdj: 0.02,
					negLog10Pval: 2,
					negLog10PvalAdj: 1.7,
					effectSize: 1.5,
				},
			],
			effectSizeLabel: "log2fc",
			skippedGenes: 0,
			elapsed: 10,
		} as any);
		vi.spyOn(DGERunner.prototype, "destroy").mockImplementation(() => {});

		const chartManager: any = {
			dsIndex: {
				cells: { dataStore: cellsDs, name: "cells" },
				rna: { dataStore: genesRna, name: "rna" },
				protein: { dataStore: genesProtein, name: "protein" },
			},
			loadColumnSetAsync: vi.fn(async (columns: string[]) => {
				for (const field of columns) {
					if (!cellsDs.columnIndex[field]) {
						const buf = new SharedArrayBuffer(nCells * 4);
						new Float32Array(buf).set([0, 1, 2]);
						cellsDs.columnIndex[field] = { buffer: buf };
					}
				}
			}),
		};

		await runDGEOnDataStore(chartManager, "cells", "protein", {
			groupColumn: "group",
			targetGroup: "target",
			referenceGroup: "ref",
		});

		const esCol = genesProtein.columnIndex.dge_effect_size.data as Float32Array;
		expect(esCol[0]).toBeCloseTo(1.5, 5);
		expect(esCol[2]).toBeNaN();
		expect(genesRna.columnIndex.dge_effect_size).toBeUndefined();
	});

	test("throws clear error when selected genes datasource is not linked", async () => {
		const nCells = 2;
		const groupBuffer = new SharedArrayBuffer(nCells);
		new Uint8Array(groupBuffer).set([0, 1]);
		const cellsDs: any = {
			size: nCells,
			name: "cells",
			columnIndex: {
				group: { datatype: "text", values: ["A", "B"], buffer: groupBuffer },
			},
			filterBuffer: new SharedArrayBuffer(nCells),
		};
		const genesRna = makeGenesDs(2);

		vi.spyOn(linkUtils, "getRowsAsColumnsLinks").mockReturnValue([
			{
				linkedDs: { name: "rna", dataStore: { columnIndex: { gene_name: { data: new Uint8Array([0, 1]) } } } },
				link: {
					name_column: "gene_name",
					subgroups: { sg_rna: { type: "sparse", name: "RNA", label: "RNA" } },
					valueToRowIndex: new Map<string, number>(),
					initPromise: Promise.resolve(),
				},
			},
		] as any);

		const chartManager: any = {
			dsIndex: {
				cells: { dataStore: cellsDs, name: "cells" },
				rna: { dataStore: genesRna, name: "rna" },
				protein: { dataStore: makeGenesDs(2), name: "protein" },
			},
			loadColumnSetAsync: vi.fn(async () => {}),
		};

		await expect(
			runDGEOnDataStore(chartManager, "cells", "protein", {
				groupColumn: "group",
				targetGroup: "A",
				referenceGroup: "B",
			}),
		).rejects.toThrow('No rows_as_columns link found from "cells" to "protein"');
	});
});
