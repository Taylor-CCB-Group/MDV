import {
    type GenomeLocation,
    type GenomeLocationValue,
    type GenomeViewMargins,
    applyViewMargins,
    getLocationFieldsFromGenome,
    locationFromFieldValues,
} from "../genomicLocationUtils";

export type IGVBrowserSingleLocation = GenomeLocation;
export type IGVBrowserLocation = GenomeLocationValue;

export type IGVBrowserViewMargins = GenomeViewMargins;

export interface IGVBaseFeature {
    chr: string;
    start: number;
    end: number;
    id: string;
    color?: string;
    name?: string;
    chr2?: string;
    pos1?: number;
    pos2?: number;
    svtype?: string;
    length?: number;
    _f?: IGVBaseFeature;//sometimes igv wraps features in an object with a _f property
}

export interface IGVStructuralVariantStyle {
    fillStyle: string;
    strokeStyle: string;
    glyph: "capped_line" | "diamond" | "triangle" | "double_bar" | "inversion" | "breakend";
}

export interface BuildBaseFeaturesResult {
    features: IGVBaseFeature[];
}

export function getStructuralVariantStyle(svtype: unknown): IGVStructuralVariantStyle {
    const normalized = typeof svtype === "string" ? svtype.trim().toUpperCase() : "";
    switch (normalized) {
        case "DEL":
            return { fillStyle: "#d94841", strokeStyle: "#922b21", glyph: "capped_line" };
        case "INS":
            return { fillStyle: "#8e44ad", strokeStyle: "#5b2c6f", glyph: "triangle" };
        case "INV":
            return { fillStyle: "#2ecc71", strokeStyle: "#1e8449", glyph: "inversion" };
        case "DUP":
            return { fillStyle: "#2980b9", strokeStyle: "#1f618d", glyph: "double_bar" };
        case "TRA":
        case "BND":
            return { fillStyle: "#e67e22", strokeStyle: "#af601a", glyph: "breakend" };
        default:
            return { fillStyle: "#777777", strokeStyle: "#555555", glyph: "diamond" };
    }
}

export function buildBaseFeatures(
    rows: Record<string, unknown>[],
    locationFields: string[],
    isSvs: boolean,
    maxInitialFeatures: number,
): BuildBaseFeaturesResult {
    if (locationFields.length < 3) {
        return { features: [] };
    }

    const features: IGVBaseFeature[] = [];
    const maxRows = Math.min(rows.length, maxInitialFeatures);

    for (let i = 0; i < maxRows; i++) {
        const row = rows[i];
        if (!row) continue;

        const chr = row[locationFields[0]];
        const startValue = Number(row[locationFields[1]]);
        const endValue = Number(row[locationFields[2]]);
        const rowId = String((row.__index__ ?? i));
        if (typeof chr !== "string" || !Number.isFinite(startValue) || !Number.isFinite(endValue)) {
            continue;
        }

        if (isSvs && locationFields.length >= 4) {
            const chr2 = row[locationFields[3]];
            const svtype = typeof row.svtype === "string" ? row.svtype : undefined;
            const normalizedSvType = typeof svtype === "string" ? svtype.trim().toUpperCase() : "";
            const style = getStructuralVariantStyle(svtype);
            
            const baseFeature = {
                id: rowId,
                chr2: typeof chr2 === "string" ? chr2 : undefined,
                pos1: startValue,
                pos2: endValue,
                svtype,
                length: Number.isFinite(Number(row.length)) ? Number(row.length) : undefined,
                name: svtype,
                color: style.fillStyle,
            };

            if (normalizedSvType === "TRA" || normalizedSvType === "BND") {
                features.push({
                    ...baseFeature,
                    id: rowId,
                    chr,
                    start: startValue,
                    end: startValue,
                });
                if (typeof chr2 === "string") {
                    features.push({
                        ...baseFeature,
                        id: rowId,
                        chr: chr2,
                        start: endValue,
                        end: endValue,
                    });
                }
                continue;
            }

            if (normalizedSvType === "INS") {
                features.push({
                    ...baseFeature,
                    chr,
                    start: startValue,
                    end: startValue,
                });
                continue;
            }

            features.push({
                ...baseFeature,
                chr,
                start: Math.min(startValue, endValue),
                end: Math.max(startValue, endValue),
            });
            continue;
        }

        features.push({
            chr,
            start: Math.min(startValue, endValue),
            end: Math.max(startValue, endValue),
            id: rowId,
        });
    }

    return { features };
}

export { applyViewMargins, getLocationFieldsFromGenome, locationFromFieldValues };


