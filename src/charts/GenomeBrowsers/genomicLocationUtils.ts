export interface GenomeLocation {
    chr: string;
    start: number;
    end: number;
}

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
        return [cols.chr1, cols.pos1, cols.pos2, cols.chr2];
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
        chr: unknown;
        start: unknown;
        end: unknown;
        chr2?: unknown;
    },
    isSvs: boolean,
): GenomeLocation | null {
    const chr = values.chr;
    const start = values.start;
    const endCandidate = values.end;
    if (typeof chr !== "string" || typeof start !== "number" || typeof endCandidate !== "number") {
        return null;
    }

    let end = endCandidate;
    if (isSvs) {
        const chr2 = values.chr2;
        if (typeof chr2 !== "string" || chr2 !== chr || Math.abs(start - endCandidate) > 30000) {
            end = start + 1000;
        }
    }

    return { chr, start, end };
}

export function highlightedIndexToLocation(
    highlightedIndex: number,
    fieldSpecs: ValueGetterSpec[],
    isSvs: boolean,
): GenomeLocation | null {
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
            chr: chrSpec.getValue(highlightedIndex),
            start: startSpec.getValue(highlightedIndex),
            end: endSpec.getValue(highlightedIndex),
            chr2: isSvs && fieldSpecs[3] ? fieldSpecs[3].getValue(highlightedIndex) : undefined,
        },
        isSvs,
    );
}
