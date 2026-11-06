import { loadColumn } from "@/dataloaders/DataLoaderUtil";
import type DataStore from "@/datastore/DataStore";
import { type AnalysisFrame, regionsScaleFrame } from "./analysisFrame";
import type { CrossPcfRequest, CrossPcfResponse } from "./crossPcfCompute";
import type { DomainSource } from "./domainSource";
import { readRegionStore } from "./regionSpatialData";
import { spatialDataRegionFrame } from "./spatialDataFrame";

/**
 * The analysis frame for rows of `region` (a key of `regions.all_regions`): the SpatialData
 * transform of the region's image when the project was converted from SpatialData, else the
 * `regions.scale` stopgap. With no region, a datastore with a single SpatialData region uses that one
 * (converted projects also carry image-only regions such as a lossy OME-TIFF, which don't count).
 */
export async function resolveAnalysisFrame(dataStore: DataStore, region?: string): Promise<AnalysisFrame> {
    const fallback = regionsScaleFrame(dataStore.regions);
    const { entry, sdata, error } = await readRegionStore(dataStore, region);
    if (error) return { ...fallback, provenance: `${fallback.provenance} (${error})` };
    if (!entry || !sdata) return fallback;
    const frame = spatialDataRegionFrame(sdata, entry);
    if ("error" in frame) return { ...fallback, provenance: `${fallback.provenance} (${frame.error})` };
    return frame;
}

export type CrossPcfOptions = {
    /** Categorical column whose categories are the cell types. Default: `regions.default_color` if
     *  categorical, else the first text column with 2–256 categories. */
    typeColumn?: string;
    /** Default: `resolveAnalysisFrame(dataStore, source.region)`. Domain and cells both use it. */
    frame?: AnalysisFrame;
    /** In the frame's unit. */
    rMax?: number;
    nBins?: number;
    /** Clip each anchor's annulus to the domain — exact near the boundary. Default true. */
    edgeCorrected?: boolean;
    /** Abort terminates the worker; the promise rejects with the signal's reason. */
    signal?: AbortSignal;
};

export type CrossPcfResult = CrossPcfResponse & {
    source: string;
    typeColumn: string;
    /** Category labels, index-aligned to the matrix axes. */
    labels: string[];
    frame: AnalysisFrame;
};

/**
 * Cross-PCF matrix of the cells inside a domain, with the domain as the observation window
 * (intraspatial `crossPCFMatrixBinned`), computed in a worker. The one entry point, whatever the
 * domain came from.
 */
export async function computeCrossPCF(
    dataStore: DataStore,
    source: DomainSource,
    opts: CrossPcfOptions = {},
): Promise<CrossPcfResult> {
    const regions = dataStore.regions;
    const typeColumn = opts.typeColumn ?? defaultTypeColumn(dataStore);
    if (!typeColumn) throw new Error("No categorical cell-type column found: pass typeColumn");
    const regionField = source.region !== undefined ? regions?.region_field : undefined;
    const frame = opts.frame ?? (await resolveAnalysisFrame(dataStore, source.region));

    const [xCol, yCol, typeCol, regionCol] = await Promise.all([
        loadColumn(dataStore.name, source.columns[0]),
        loadColumn(dataStore.name, source.columns[1]),
        loadColumn(dataStore.name, typeColumn),
        regionField ? loadColumn(dataStore.name, regionField) : undefined,
    ]);
    const labels = typeCol.values;
    if (!labels) throw new Error(`${typeColumn} is not a categorical column`);

    const request: CrossPcfRequest = {
        x: { buffer: xCol.data.buffer, datatype: xCol.datatype },
        y: { buffer: yCol.data.buffer, datatype: yCol.datatype },
        type: { buffer: typeCol.data.buffer, datatype: typeCol.datatype, nTypes: labels.length },
        frame,
        domain: await source.toDomain(frame),
        rMax: opts.rMax ?? 100,
        nBins: opts.nBins ?? 20,
        edgeCorrected: opts.edgeCorrected ?? true,
    };
    if (regionCol && source.region !== undefined) {
        const code = regionCol.values?.indexOf(source.region) ?? -1;
        if (code < 0) throw new Error(`Region ${source.region} not found in ${regionField}`);
        request.region = { buffer: regionCol.data.buffer, datatype: regionCol.datatype, code };
    }

    const worker = new Worker(new URL("./crossPcfWorker.ts", import.meta.url), { type: "module" });
    try {
        const result = await new Promise<CrossPcfResponse>((resolve, reject) => {
            opts.signal?.throwIfAborted();
            opts.signal?.addEventListener("abort", () => reject(opts.signal?.reason), { once: true });
            worker.onmessage = (e: MessageEvent<{ result?: CrossPcfResponse; error?: string }>) =>
                e.data.result ? resolve(e.data.result) : reject(new Error(e.data.error));
            worker.onerror = (e) => reject(new Error(e.message));
            worker.postMessage(request);
        });
        return { ...result, source: source.label, typeColumn, labels, frame };
    } finally {
        worker.terminate();
    }
}

function defaultTypeColumn(dataStore: DataStore): string | undefined {
    const isCategorical = (c: { datatype: string; values?: unknown[] } | undefined) =>
        !!c &&
        (c.datatype === "text" || c.datatype === "text16") &&
        (c.values?.length ?? 0) >= 2 &&
        (c.values?.length ?? 0) <= 256;
    const preferred = dataStore.regions?.default_color;
    if (preferred && isCategorical(dataStore.columnIndex[preferred])) return preferred;
    return dataStore.columns.find((c) => isCategorical(c) && c.field !== dataStore.regions?.region_field)?.field;
}

/** Debug output: the r-averaged g for every non-empty pair. */
export function logCrossPCF(res: CrossPcfResult) {
    const { nTypes: K, nBins: B, g, counts, labels } = res;
    const present = [...Array(K).keys()].filter((t) => counts[t] > 0);
    const meanG: Record<string, Record<string, string>> = {};
    for (const a of present) {
        const row: Record<string, string> = {};
        for (const b of present) {
            let s = 0;
            for (let k = 0; k < B; k++) s += g[a * K * B + b * B + k];
            row[labels[b]] = (s / B).toFixed(2);
        }
        meanG[`${labels[a]} (n=${counts[a]})`] = row;
    }
    const inside = res.candidates - (res.outside ?? 0);
    console.log(
        `[cross-PCF] ${res.source}: ${inside} cells inside (of ${res.candidates} in bbox), ` +
            `${present.length}/${K} ${res.typeColumn} types, r ∈ [0, ${res.r[B - 1] + res.r[0]}) ${res.frame.unit}, ` +
            `${B} bins, ${res.ms.toFixed(0)} ms; frame: ${res.frame.provenance}`,
    );
    console.log("[cross-PCF] mean g over r (row = anchor type, column = neighbour type)");
    console.table(meanG);
    return meanG;
}
