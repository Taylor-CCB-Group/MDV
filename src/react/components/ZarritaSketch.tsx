import * as zarr from "zarrita";
import { ZipFileStore } from "@zarrita/storage";
import { useEffect, useState } from "react";
import JsonView from "react18-json-view";

function makeXeniumCellId(cellIdPrefix: number, datasetSuffix: number) {
    // Step 1: Convert to hex and pad to 8 characters
    const hex = cellIdPrefix.toString(16).padStart(8, '0');

    // Step 2: Shift hex characters to the custom range a–p
    const shiftMap = {
        '0': 'a', '1': 'b', '2': 'c', '3': 'd',
        '4': 'e', '5': 'f', '6': 'g', '7': 'h',
        '8': 'i', '9': 'j', 'a': 'k', 'b': 'l',
        'c': 'm', 'd': 'n', 'e': 'o', 'f': 'p'
    } as const;

    //@ts-expect-error - could have a better type for shiftMap, or verify that hex only contains valid hex characters
    const shifted = hex.split('').map(char => shiftMap[char]).join('');

    // Step 3: Append dash and suffix
    return `${shifted}-${datasetSuffix}`;
}

function parseXeniumCellId(cellId: string) {
    // Step 1: Split into shifted prefix and suffix
    const [shifted, suffixStr] = cellId.split('-');

    // Step 2: Reverse the custom shift mapping
    const reverseMap = {
        'a': '0', 'b': '1', 'c': '2', 'd': '3',
        'e': '4', 'f': '5', 'g': '6', 'h': '7',
        'i': '8', 'j': '9', 'k': 'a', 'l': 'b',
        'm': 'c', 'n': 'd', 'o': 'e', 'p': 'f'
    };

    // @ts-expect-error - could have a better type for reverseMap, or verify that shifted only contains valid characters
    const hex = shifted.split('').map(char => reverseMap[char]).join('');

    // Step 3: Convert hex string to integer
    const cell_id_prefix = Number.parseInt(hex, 16);

    // Step 4: Parse suffix
    const dataset_suffix = Number.parseInt(suffixStr, 10);

    return { cell_id_prefix, dataset_suffix };
}

function reconstructCellIds(cellIdArray: zarr.TypedArray<zarr.NumberDataType>) {
    const n = cellIdArray.length / 2;
    const ids = new Array(n);
    for (let i = 0; i < n; i++) {
        const prefix = cellIdArray[i * 2];
        const suffix = cellIdArray[i * 2 + 1];
        ids[i] = makeXeniumCellId(prefix, suffix);
    }
    return ids;
}

const url = "http://localhost:8080";
async function test() {
    // these xenium stores don't have the 'consolidated' metadata 
    // - but if we're making our own spec, we likely will.
    // const cellsStore = await zarr.tryWithConsolidated(new zarr.FetchStore(`${url}/cells.zarr`));
    const cellsStore = await zarr.tryWithConsolidated(ZipFileStore.fromUrl(`${url}/cells.zarr.zip`)) as zarr.Readable;
    const analysisStore = await zarr.tryWithConsolidated(ZipFileStore.fromUrl(`${url}/analysis.zarr.zip`)) as zarr.Readable;
    // const transcriptsStore = await zarr.tryWithConsolidated(ZipFileStore.fromUrl(`${url}/transcripts.zarr.zip`)) as zarr.Readable;

    const cellsRoot = await zarr.open(cellsStore);
    const anRoot = await zarr.open(analysisStore);

    const id = await zarr.open(cellsRoot.resolve("cell_id"), { kind: "array" });
    if (!id.is("number")) throw new Error("Expected cell_id to be an array of numbers");

    const summary = await zarr.open(cellsRoot.resolve("cell_summary"), { kind: "array"});
    if (!summary.is("number")) throw new Error("Expected cell_summary to be an array of numbers");
    const ids = await zarr.get(id);
    // can we find a specific cell `kalnagga-1`?
    // it should have location 665.253552.41 3552.41, area 87.3, Cluster 8...
    // let's find the indices & indptrs for cluster 8 (from Graph-Based Clustering, ie group 0)...
    const clusterIndices = await zarr.open(anRoot.resolve("cell_groups/0/indices"), { kind: "array" });
    if (!clusterIndices.is("number")) throw new Error("Expected cluster indices to be an array of numbers");
    const clusterIndptr = await zarr.open(anRoot.resolve("cell_groups/0/indptr"), { kind: "array" });
    if (!clusterIndptr.is("number")) throw new Error("Expected cluster indptr to be an array of numbers");
    // I should be able to find the indptr (7, zero-based), to know where to start looking in indices for our cluster... I'll also want the next one so I know the region to slice
    const indPtrs = await zarr.get(clusterIndptr); //could have sliced here, but meh... it's small
    
    const c8Slice = zarr.slice(indPtrs.data[7], indPtrs.data[8]);
    const c8Indices = await zarr.get(clusterIndices, [c8Slice]); //this gives 2966 indices, as expected

    //... so I should be able to use c8Indices to... profit???

    // this is heavy, we'd likely want to do things differently
    const allCellIdStrings = reconstructCellIds(ids.data);
    //still trying to find this one cell I decided I was interested in...
    //we also know the area is ~87.3µm^2 - `/cell_summary` should have `[x, y, cell_area, ...]`...
    const allCellAreas = await zarr.get(summary, [null, 2]);
    //! oops... if we map over a Uint32Array, we get another Uint32Array... 
    //but we should have a Float64Array here
    // const c8Areas = new Float64Array(c8Indices.data).map(i => allCellAreas.data[i]);
    const c8Areas = Float64Array.from(c8Indices.data, i => allCellAreas.data[i]);
    const queriedCellSize = 87.3;
    // find the indices of cells with areas similar to what we're looking for...
    const closeIndices = c8Areas.map((v, i) => Math.abs(queriedCellSize - v) <= 0.1 ? c8Indices.data[i] : -1).filter(i => i !== -1);
    // now let's get the centroids for these cells...
    // const queryX = 665.25, queryY = 3552.41; //for ref, this is the one we are hoping to find
    const allCellCentroids = await zarr.get(summary, [null, zarr.slice(0, 2)]);
    // or is there a transform matrix involved?
    const transform = await zarr.open(cellsRoot.resolve("/masks/homogeneous_transform"), { kind: "array" });
    // looks like a scale matrix, with x&y 4.7058820724487305
    // const transformData = await zarr.get(transform, [null, null]);


    const numVerts = await zarr.open(cellsRoot.resolve("polygon_num_vertices"), { kind: "array" });
    if (!numVerts.is("number")) throw new Error("Expected polygon_num_vertices to be an array of numbers");
    const numVertsArr = await zarr.get(
        numVerts,
        [0, null] // shape: [2, N] do we have a both nucleus and cell boundary in here maybe?
    );
    const numVertsSet = new Set(Array.from(numVertsArr.data));
    function computeOffsets(numVerts: any) {
        if (numVerts.length === 0) {
            return new Uint32Array(0);
        }
        const offsets = new Uint32Array(numVerts.length);
        let sum = 0;
        for (let i = 0; i < numVerts.length; i++) {
            offsets[i] = sum;
            sum += numVerts[i];
        }
        return offsets;
    }

    const offsets = computeOffsets(numVertsArr.data);
    
    // broken...
    async function loadPolygonVertices(cellIndex: number) {
        const start = offsets[cellIndex];

        const count = numVertsArr.data[cellIndex];
        const polyVertsArr = await zarr.open(cellsRoot.resolve("polygon_vertices"), { kind: "array" });
        if (!polyVertsArr.is("number")) throw new Error("Expected polygon_vertices to be an array of numbers");
        try {
            const polyVerts = await zarr.get(
                polyVertsArr,
                [zarr.slice(0, 2), zarr.slice(start, start + count)]
            );
            return {
                x: polyVerts.data.subarray(0, count),
                y: polyVerts.data.subarray(count)
            };
        } catch (error) {
            return error;
        }
    }


    const queriedCells = await Promise.all([...closeIndices].map(async i => {
        const start = offsets[i];
        const count = numVertsArr.data[i]; //always 13? Is that right? Unlucky for some...
        
        const polyVerts = await loadPolygonVertices(i); //still broken...
        return {
            id: allCellIdStrings[i],
            x: allCellCentroids.data[i * 2 + 0],
            y: allCellCentroids.data[i * 2 + 1],
            cellArea: allCellAreas.data[i],
            start, count,
            polyVerts
        }
    }));

    // numVertsArr.data is Uint32Array of shape [nCells]
    // const summary0 = (await zarr.get(summary, [zarr.slice(1), null]));
    // // how many vertices does this cell have?
    // const polyNumVerts = await zarr.open(cellsRoot.resolve("polygon_num_vertices"));
    const polyVerts = await zarr.open(cellsRoot.resolve("polygon_vertices"));
    // const numVerts0 = (await zarr.get(polyNumVerts, [zarr.slice(1), 0])).data[0];
    // // a whole lot of 0?
    // const verts = (await zarr.get(polyVerts, [1, zarr.slice(0, numVerts0)]));

    return {
        queriedCells, offsets, numVertsSet, polyVerts, numVertsArr
    }
}


async function transcriptTest() {
    const transcriptsStore = await zarr.tryWithConsolidated(ZipFileStore.fromUrl(`${url}/transcripts.zarr.zip`)) as zarr.Readable;
    // let's try to get the locations of all transcripts in the most zoomed out pyramid level, with a quality threshold
    const qualityThreshold = 10;
    const transcriptsRoot = await zarr.open(transcriptsStore);
    // todo get the grids & check the shape to find the most zoomed out level...
    // for now, cheating & hardcoding the path... in a spatial view, we'll need to resolve pyramid levels & locations properly
    const grids = await zarr.open(transcriptsRoot.resolve("grids/6/0,0"), { kind: "group" });
    const location = await zarr.open(grids.resolve("location"), { kind: "array" });
    if (!location.is("number")) throw new Error("Expected location to be an array of numbers");
    const quality = await zarr.open(grids.resolve("quality_score"), { kind: "array" });
    if (!quality.is("float32")) throw new Error("Expected quality_scores to be float32");
    // let's do that quality filter - find indices where quality is above 30
    const qualityData = await zarr.get(quality);
    const highQualityIndices = Uint32Array.from(
        Array.from(qualityData.data).map((q, i) => q <= qualityThreshold ? i : -1).filter(i => i !== -1)
    );
    // now let's get the locations for these indices
    const locationData = await zarr.get(location);
    const locations = Array.from(highQualityIndices).map(i => {
        const x = locationData.data[i * 3];
        const y = locationData.data[i * 3 + 1];
        const z = locationData.data[i * 3 + 2];
        const q = qualityData.data[i];
        return { x, y, z, q };
    });
    
    return { locations };
}

export default function ZarritaSketch() {
    const url = "http://localhost:8080/";
    const [view, setView] = useState<any>();
    useEffect(() => {
        transcriptTest().then(setView);
    }, []); // url should be defined in the parent component


    return <JsonView src={view} />;
}