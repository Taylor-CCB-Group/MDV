/**
 * Bridges DGERunner with MDV's live ChartManager/DataStore infrastructure.
 *
 * Resolves gene fields from rows_as_columns links, loads columns in batches,
 * runs the DGE pipeline, and writes result columns back to the genes DataStore.
 */

import type { DGEConfig, DGERunResult, EffectSizeLabel } from "./DGEDimension";
import { DGERunner } from "./DGEDimension";
import { detectDataIsLog1p, clampEffectSize } from "./dgeStats";
import { getRowsAsColumnsLinks, getFieldName } from "../links/link_utils";

type ChartManager = {
	dsIndex: Record<string, { dataStore: any; name: string } | undefined>;
	dataSources: { name: string; dataStore: any }[];
	loadColumnSetAsync(columns: string[], dataSource: string): Promise<void>;
	addChart(dataSource: string, config: Record<string, unknown>, notify?: boolean): Promise<void>;
};

export interface DGEIntegrationConfig {
	groupColumn: string;
	targetGroup: string;
	referenceGroup: string | "rest";
	batchSize?: number;
}

export interface DGEIntegrationResult extends DGERunResult {
	cellsDsName: string;
	genesDsName: string;
}

/**
 * Run DGE analysis using live MDV DataStore infrastructure.
 *
 * 1. Resolves gene fields from the rows_as_columns link
 * 2. Detects whether expression data is log1p or z-scored
 * 3. Runs the DGE pipeline via DGERunner
 * 4. Writes result columns to the genes DataStore
 */
export async function runDGEOnDataStore(
	chartManager: ChartManager,
	cellsDsName: string,
	genesDsName: string,
	config: DGEIntegrationConfig,
	onProgress?: (done: number, total: number) => void,
): Promise<DGEIntegrationResult> {
	const cellsDs = chartManager.dsIndex[cellsDsName]?.dataStore;
	const genesDs = chartManager.dsIndex[genesDsName]?.dataStore;
	if (!cellsDs) throw new Error(`DataStore "${cellsDsName}" not found`);
	if (!genesDs) throw new Error(`DataStore "${genesDsName}" not found`);

	// Load group column
	await chartManager.loadColumnSetAsync([config.groupColumn], cellsDsName);
	const groupCol = cellsDs.columnIndex[config.groupColumn];
	if (!groupCol) throw new Error(`Group column "${config.groupColumn}" not found`);
	if (!groupCol.values) throw new Error(`Group column "${config.groupColumn}" is not categorical`);

	// Normalize group assignments into a Uint8Array SharedArrayBuffer.
	// Text columns use Uint8Array, text16 columns use Uint16Array — we need
	// a consistent Uint8Array for the DGE worker regardless of source type.
	const nCells = cellsDs.size as number;
	const groupSAB = new SharedArrayBuffer(nCells);
	const groupArr = new Uint8Array(groupSAB);
	if (groupCol.datatype === "text16") {
		const src = new Uint16Array(groupCol.buffer);
		for (let i = 0; i < nCells; i++) groupArr[i] = src[i];
	} else {
		const src = new Uint8Array(groupCol.buffer);
		groupArr.set(src);
	}

	// Resolve gene fields from the rows_as_columns link
	const racLinks = getRowsAsColumnsLinks(cellsDs);
	if (racLinks.length === 0) {
		throw new Error(`No rows_as_columns link found for "${cellsDsName}"`);
	}
	const { linkedDs, link } = racLinks[0];
	await link.initPromise;

	if (!link.valueToRowIndex) {
		throw new Error("Link valueToRowIndex not initialized");
	}

	// Build gene field names and corresponding gene names from the linked DataStore
	const nameCol = linkedDs.dataStore.columnIndex[link.name_column];
	if (!nameCol?.data) {
		throw new Error(`Name column "${link.name_column}" not loaded in genes DataStore`);
	}

	const sgKey = Object.keys(link.subgroups)[0];
	const geneNames: string[] = [];
	const geneFields: string[] = [];
	const geneRowIndices: number[] = [];

	for (const [geneName, rowIndex] of link.valueToRowIndex) {
		geneNames.push(geneName);
		geneFields.push(getFieldName(sgKey, geneName, rowIndex));
		geneRowIndices.push(rowIndex);
	}

	// ── Diagnostics: group column ──
	const groupValues: string[] = groupCol.values;
	const targetIdx = groupValues.indexOf(config.targetGroup);
	const refIdx = config.referenceGroup === "rest" ? -1 : groupValues.indexOf(config.referenceGroup);
	const groupCounts = new Map<number, number>();
	for (let i = 0; i < nCells; i++) {
		groupCounts.set(groupArr[i], (groupCounts.get(groupArr[i]) ?? 0) + 1);
	}
	const filterArr = new Uint8Array(cellsDs.filterBuffer);
	let filteredCount = 0;
	for (let i = 0; i < nCells; i++) if (filterArr[i] !== 0) filteredCount++;
	console.log("[DGE diag] Group column:", config.groupColumn, "datatype:", groupCol.datatype);
	console.log("[DGE diag] Group values:", groupValues);
	console.log("[DGE diag] Target:", config.targetGroup, `(idx=${targetIdx})`, "Reference:", config.referenceGroup, `(idx=${refIdx})`);
	console.log("[DGE diag] Group distribution:", Object.fromEntries([...groupCounts].map(([k, v]) => [groupValues[k] ?? `?${k}`, v])));
	console.log("[DGE diag] Total cells:", nCells, "Filtered out:", filteredCount, "Active:", nCells - filteredCount);
	console.log("[DGE diag] Total genes:", geneFields.length);

	// Determine whether expression data is log1p-normalized or z-scored.
	// Prefer the subgroup storage type from datasources.json metadata:
	//   sparse → data preserves zero structure → log1p (z-scoring densifies)
	//   dense / missing → ambiguous, fall back to probing gene values
	const sgType = link.subgroups[sgKey]?.type as string | undefined;
	let dataIsLog1p = true;

	if (sgType === "sparse") {
		dataIsLog1p = true;
		console.log("[DGE diag] sgtype=sparse -> dataIsLog1p=true (from metadata)");
	} else {
		console.log(`[DGE diag] sgtype=${sgType ?? "undefined"} -> probing genes for detection`);
		const MAX_PROBE = Math.min(20, geneFields.length);
		const probeFields = geneFields.slice(0, MAX_PROBE);
		await chartManager.loadColumnSetAsync(probeFields, cellsDsName);

		for (let gi = 0; gi < probeFields.length; gi++) {
			const col = cellsDs.columnIndex[probeFields[gi]];
			const buf = col?.buffer as SharedArrayBuffer | null;
			if (!buf) continue;
			const arr = new Float32Array(buf);
			const detection = detectDataIsLog1p(arr);
			if (detection !== null) {
				dataIsLog1p = detection;
				console.log(`[DGE diag] Probed gene #${gi} "${geneNames[gi]}" -> dataIsLog1p=${dataIsLog1p}`);
				break;
			}
			if (gi === probeFields.length - 1) {
				console.warn(`[DGE diag] All ${MAX_PROBE} probed genes were all-NaN; defaulting to dataIsLog1p=true`);
			}
		}
	}

	const dgeConfig: DGEConfig = {
		groupColumn: config.groupColumn,
		targetGroup: config.targetGroup,
		referenceGroup: config.referenceGroup,
		// Adjust this parameter to control the DGE batch size (number of genes per request).
		// Higher values (e.g. 500-1000) are faster but use more memory per batch.
		batchSize: config.batchSize ?? 2000,
	};

	const runner = new DGERunner();
	try {
		const result = await runner.run(
			dgeConfig,
			cellsDs.filterBuffer as SharedArrayBuffer,
			nCells,
			groupSAB,
			groupValues,
			geneFields,
			geneNames,
			async (fields: string[]) => {
				await chartManager.loadColumnSetAsync(fields, cellsDsName);
			},
			(field: string) => {
				const col = cellsDs.columnIndex[field];
				return (col?.buffer as SharedArrayBuffer) ?? null;
			},
			onProgress,
			false,
			dataIsLog1p,
		);

		writeResultsToGenesDS(genesDs, result, link.valueToRowIndex);

		return {
			...result,
			cellsDsName,
			genesDsName,
		};
	} finally {
		runner.destroy();
	}
}

/**
 * Write DGE result columns into the genes DataStore, mapping results back
 * to the correct row order using the gene name → row index map.
 */
function writeResultsToGenesDS(
	genesDs: any,
	result: DGERunResult,
	valueToRowIndex: Map<string, number>,
): void {
	const n = genesDs.size as number;
	const esLabel = result.effectSizeLabel;
	const esColumnName = esLabel === "log2fc" ? "dge_log2fc" : "dge_mean_diff";

	const columns: { field: string; name: string; data: Float32Array }[] = [
		{ field: esColumnName, name: esLabel === "log2fc" ? "Log2 FC" : "Mean Diff", data: new Float32Array(n).fill(NaN) },
		{ field: "dge_effect_size", name: "Effect Size", data: new Float32Array(n).fill(NaN) },
		{ field: "dge_pval", name: "P-value", data: new Float32Array(n).fill(NaN) },
		{ field: "dge_pval_adj", name: "Adj. P-value", data: new Float32Array(n).fill(NaN) },
		{ field: "dge_neg_log10_pval", name: "-log10(p)", data: new Float32Array(n).fill(NaN) },
		{ field: "dge_neg_log10_pval_adj", name: "-log10(padj)", data: new Float32Array(n).fill(NaN) },
	];

	for (const gene of result.results) {
		const rowIdx = valueToRowIndex.get(gene.gene);
		if (rowIdx === undefined) continue;

		const clampedES = clampEffectSize(gene.effectSize);
		columns[0].data[rowIdx] = clampedES;
		columns[1].data[rowIdx] = clampedES;
		columns[2].data[rowIdx] = gene.pval;
		columns[3].data[rowIdx] = gene.pvalAdj;
		columns[4].data[rowIdx] = gene.negLog10Pval;
		columns[5].data[rowIdx] = gene.negLog10PvalAdj;
	}

	const addedFields: string[] = [];
	for (const col of columns) {
		const buffer = new SharedArrayBuffer(col.data.byteLength);
		new Float32Array(buffer).set(col.data);

		if (genesDs.columnIndex[col.field]) {
			genesDs.setColumnData(col.field, buffer);
		} else {
			genesDs.addColumn(
				{ name: col.name, field: col.field, datatype: "double" },
				buffer,
				true,
			);
		}
		addedFields.push(col.field);
	}

	genesDs.dataChanged(addedFields, true);
}

/**
 * Find which datasource pairs have a rows_as_columns link,
 * returning the cells DS name and linked genes DS name.
 */
export function findDGECapableDatasources(
	chartManager: ChartManager,
): { cellsDsName: string; genesDsName: string }[] {
	const results: { cellsDsName: string; genesDsName: string }[] = [];
	for (const ds of chartManager.dataSources) {
		const dataStore = ds.dataStore;
		if (!dataStore.links) continue;
		for (const linkedDsName of Object.keys(dataStore.links)) {
			const linkSet = dataStore.links[linkedDsName];
			if (linkSet.rows_as_columns) {
				results.push({
					cellsDsName: ds.name,
					genesDsName: linkedDsName,
				});
			}
		}
	}
	return results;
}
