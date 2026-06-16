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
        const chr = typeof row[map.chr] === "string" ? (row[map.chr] as string).trim() : "";
        const start = Number(row[map.start]);
        const end = Number(row[map.end]);
        if (!chr || !Number.isFinite(start) || !Number.isFinite(end)) {
            return null;
        }
        return [{
            chr,
            start,
            end,
        }];
    }
    else if (genomicInfo.type === "sv") {
        const map = genomicInfo.columns ;
        const svtype = typeof row[map.svtype] === "string" ? (row[map.svtype] as string).trim() : "";
        const chr1 = typeof row[map.chr1] === "string" ? (row[map.chr1] as string).trim() : "";
        const chr2 = typeof row[map.chr2] === "string" ? (row[map.chr2] as string).trim() : "";
        const pos1 = Number(row[map.pos1]);
        const pos2 = Number(row[map.pos2]);
        if (!svtype || !chr1 || !Number.isFinite(pos1)) {
            return null;
        }
        if ((svtype === "BND" || svtype === "TRA") || pos2-pos1 > 1000000){
            if (!chr2 || !Number.isFinite(pos2)) {
                return null;
            }
            return [
                {
                    chr: chr1,
                    start: pos1,
                    end: pos1,
                },
                {
                    chr: chr2,
                    start: pos2,
                    end: pos2,
                }
            ]
        }
        else {
            if (!Number.isFinite(pos2)) {
                return null;
            }
            return [{
                chr: chr1,
                start: pos1,
                end: pos2,
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
                for (let i=0 ; i < locationFields.length; i++) {
                   if (locationColumns[i] !== undefined) {
                       columns[names[i]] = locationColumns[i];
                   }
                }
                return columns;
            }
            return {};
        },
        [areColumnsLoaded, genome, locationFields, locationColumns],
    );
    return {
        genomicInfo:dataStore.genome,
        genomicColumns,
        areColumnsLoaded,
    };
}


