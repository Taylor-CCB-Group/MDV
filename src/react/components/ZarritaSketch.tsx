import * as zarr from "zarrita";
// import { ZipFileStore } from "@zarrita/storage";
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

async function testMerfishSpatialData() {
    // this only worked after I renamed `metadata` to `.metadata`???
    const url = "http://localhost:8081/data.zarr/";
    const store = new zarr.FetchStore(url);
    const merfish = await zarr.tryWithConsolidated(store); //returns `any` according to lsp???
    // const url = "http://localhost:8081/merfish.zip";
    // const store = ZipFileStore.fromUrl(url);
    // const merfish = await zarr.tryWithConsolidated(store) as zarr.AsyncReadable; //missing key /zarr.json
    const root = await zarr.open(merfish, { kind: "group" });
    
    return { root };
}

async function testXeniumS3() {
    const url = "https://s3.embl.de/spatialdata/spatialdata-sandbox/xenium_rep2_io.zarr/";
    const store = new zarr.FetchStore(url, {
        overrides: {
            mode: "cors"
        }
    });
    const xenium = await zarr.tryWithConsolidated(store);
    return { xenium };
    // const root = await zarr.open(xenium);
}

export default function ZarritaSketch() {
    const url = "http://localhost:8080/";
    const [view, setView] = useState<any>();
    useEffect(() => {
        // transcriptTest().then(setView);
        testMerfishSpatialData().then(setView);
    }, []); // url should be defined in the parent component


    return <JsonView src={view} />;
}