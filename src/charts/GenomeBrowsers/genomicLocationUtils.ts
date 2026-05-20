export interface GenomeLocation {
    chr: string;
    start: number;
    end: number;
}

export type GenomeLocationPair = [GenomeLocation] | [GenomeLocation, GenomeLocation];
export type GenomeLocationValue = GenomeLocation | GenomeLocationPair;

export interface GenomeViewMargins {
    type: "fixed_length" | "absolute" | "percentage";
    value: number;
}

type ValueGetterSpec = {
    getValue: (rowIndex: number) => unknown;
};

export function getLocationFieldsFromGenome(genome: any): string[] | null {
    if (genome?.genomic_location?.columns) {
        const cols = genome.genomic_location.columns;
        return [cols.chr, cols.start, cols.end];
    }
    if (genome?.svs?.sv_columns) {
        const cols = genome.svs.sv_columns;
        return [cols.chr1, cols.pos1, cols.pos2, cols.chr2, cols.svtype,cols.length];
    }
    return null;
}

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
        case "absolute": {
            start -= vm.value;
            end += vm.value;
            break;
        }
        case "percentage": {
            const length = Math.max(1, end - start);
            const margin = Math.round((length * vm.value) / 100);
            start -= margin;
            end += margin;
            break;
        }
        case "fixed_length": {
            start -= vm.value;
            end += vm.value;
            break;
        }
    }
    return { chr, start: Math.max(0, start), end: Math.max(1, end) };
}

export function locationFromFieldValues(
    values: {
        chr: string;
        start: number;
        end: number;
        chr2?: string;
    },
    isSvs: boolean,
): GenomeLocationValue | null {
    const chr = values.chr;
    const start = values.start;
    const end = values.end;
  
    if (isSvs) {
    
    const chr2 = values.chr2;
       //need two locations
       if (chr2 &&chr2 !== chr) {
            const loc1: GenomeLocation = {
                chr,
                start,
                end: start
            };
            const loc2: GenomeLocation = {
                chr: chr2,
                start: end,
                end: end,
            };
            return [loc1, loc2];
        }
    }
    return { chr, start, end };
}

export function highlightedIndexToLocation(
    highlightedIndex: number,
    fieldSpecs: ValueGetterSpec[],
    isSvs: boolean,
): GenomeLocationValue | null {
    if (highlightedIndex < 0 || fieldSpecs.length < 3) {
        return null;
    }
    const chrSpec = fieldSpecs[0];
    const startSpec = fieldSpecs[1];
    const endSpec = fieldSpecs[2];
    if (!chrSpec || !startSpec || !endSpec) {
        return null;
    }
    return locationFromFieldValues(
        {
            chr: chrSpec.getValue(highlightedIndex) as string,
            start: startSpec.getValue(highlightedIndex) as number,
            end: endSpec.getValue(highlightedIndex) as number,
            chr2: isSvs && fieldSpecs[3] ? fieldSpecs[3].getValue(highlightedIndex) as string : undefined,
        },
        isSvs,
    );
}
