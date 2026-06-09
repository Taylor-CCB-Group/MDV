import { useMemo } from "react";
import { useDataStore } from "@/react/context";
import { useFieldSpecs } from "@/react/hooks";
import { DataColumn } from "../charts";


export interface GenomeLocation {
    chr: string;
    start: number;
    end: number;
}

export interface GenomeViewMargins {
    type: "fixed_length" | "absolute" | "percentage";
    value: number;
}

type GenomeMetadataBase = {
    assembly: string;
    ucsc_proxy_url?: string;
    chromosomes?: Record<string, number>;
};

type IntervalGenomeMetadata = GenomeMetadataBase & {
  type: "interval";
    columns: {
        chr: string;
        start: string;
        end: string;
    };
};

type SvGenomeMetadata = GenomeMetadataBase & {
    type: "sv";
    columns: {
        chr1: string;
        pos1: string;
        pos2: string;
        chr2: string;
        svtype: string;
        length: string;
    };
};

export type GenomeMetadata = IntervalGenomeMetadata | SvGenomeMetadata;


export function applyViewMargins(
    location: GenomeLocation,
    vm: GenomeViewMargins,
): GenomeLocation {
    const chr = location.chr;
    let start = location.start;
    let end = location.end;
    if (start > end) {
        const tmp = start;
        start = end;
        end = tmp;
    }
    switch (vm.type) {
        //add a fixed number of bases on either side
        case "absolute": {
            start -= vm.value;
            end += vm.value;
            break;
        }
        //add a percentage of the feature length on either side
        case "percentage": {
            const length = Math.max(1, end - start);
            const margin = Math.round((length * vm.value) / 100);
            start -= margin;
            end += margin;
            break;
        }
        //always the same length from the middle of the feature
        case "fixed_length": {
            const center = Math.round((start + end) / 2);   
            start = center - vm.value;
            end = center + vm.value;
            break;
        }
    }
    return { chr, start: Math.max(0, start), end: Math.max(start, end) };
}

export function locationFromFieldValues(
    dataStore:any,
    index:number,
    genomicInfo: GenomeMetadata,
): GenomeLocation[] | null {
  
    const cols = Object.values(genomicInfo.columns);
    const row = dataStore.getRowAsObject(index,cols);
    if (!row) {
        return null;
    } 
    if (genomicInfo.type === "interval") {
        const map = genomicInfo.columns ;
        return [{
            chr: row[map.chr] as string,
            start: Number(row[map.start]),
            end: Number(row[map.end]),
        }];
    }
    else if (genomicInfo.type === "sv") {
        const map = genomicInfo.columns ;
        const svtype = row[map.svtype] as string;
        const pos1 = Number(row[map.pos1]);
        const pos2 = Number(row[map.pos2])
        if ((svtype === "BND" || svtype === "TRA") || pos2-pos1 > 1000000){
            return [
                {
                    chr: row[map.chr1] as string,
                    start: Number(row[map.pos1]),
                    end: Number(row[map.pos1]),
                },
                {
                    chr: row[map.chr2] as string,
                    start: Number(row[map.pos2]),
                    end: Number(row[map.pos2]),
                }
            ]
        } 
        else {
            return [{
                chr: row[map.chr1] as string,
                start: Number(row[map.pos1]),
                end: Number(row[map.pos2]),
             }]
        }
    }
    return null;
}

export function useGenomicInfo() {
    const dataStore = useDataStore();
    const genome = dataStore.genome as GenomeMetadata;
    if (!genome) {
        throw new Error("No genome metadata found in data store");
    }
    const locationFields = useMemo(() => (Object.values(genome.columns)), [genome]);
    const locationColumns = useFieldSpecs(locationFields);
    const areColumnsLoaded =
        locationColumns.length === locationFields.length;
    const genomicColumns= useMemo(
        () => {
            if (areColumnsLoaded) {
                const names = Object.keys(genome.columns);
                const columns: Record<string,DataColumn<any>> = {};
                for (let i=0 ; i <= locationFields.length; i++) {
                   columns[names[i]] = locationColumns[i];
                }
                return columns;
            }
            return {};
        },
        [areColumnsLoaded],
    );
    return {
        genomicInfo:dataStore.genome,
        genomicColumns,
        areColumnsLoaded,
    };
}


