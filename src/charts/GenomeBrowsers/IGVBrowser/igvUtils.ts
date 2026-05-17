import {
    type GenomeLocation,
    type GenomeViewMargins,
    applyViewMargins,
    getLocationFieldsFromGenome,
    locationFromFieldValues,
} from "../genomicLocationUtils";

export type IGVBrowserLocation = GenomeLocation;

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
    rows: any[],
    columns: string[],
    isSvs: boolean,
    maxInitialFeatures: number,
) {
    const sampled = rows.length > maxInitialFeatures;
    const selectedRows = sampled ? rows.slice(0, maxInitialFeatures) : rows;

    const features: IGVBaseFeature[] = [];
    selectedRows.forEach((row, index) => {
        const chr = row?.[columns[0]];
        const startValue = Number(row?.[columns[1]]);
        const endValue = Number(row?.[columns[2]]);

        if (typeof chr !== "string" || !Number.isFinite(startValue) || !Number.isFinite(endValue)) {
            return;
        }

        const start = Math.min(startValue, endValue);
        const end = Math.max(startValue, endValue);

        if (isSvs && columns.length >= 4) {
            const chr2 = row?.[columns[3]];
            const svtype = typeof row?.svtype === "string" ? row.svtype : undefined;
            const style = getStructuralVariantStyle(svtype);
            const base = {
                svtype,
                length: Number.isFinite(Number(row?.length)) ? Number(row.length) : undefined,
                name: svtype,
                color: style.fillStyle,
                chr2: typeof chr2 === "string" ? chr2 : undefined,
                pos1: startValue,
                pos2: endValue,
                id: String(row?.__index__ ?? index),
            };
            if (svtype === "TRA" || svtype === "BND") {
                // Emit two features: one for each breakend
                features.push({
                    chr,
                    start: startValue,
                    end: startValue + 1,
                    ...base,
                    id: String(row?.__index__ ?? index) + ":1",
                } as IGVBaseFeature);
                if (typeof chr2 === "string" && Number.isFinite(endValue)) {
                    features.push({
                        chr: chr2,
                        start: endValue,
                        end: endValue + 1,
                        ...base,
                        id: String(row?.__index__ ?? index) + ":2",
                    } as IGVBaseFeature);
                }
            } else {
                features.push({
                    chr,
                    start,
                    end,
                    ...base,
                } as IGVBaseFeature);
            }
        } else {
            features.push({
                chr,
                start,
                end,
                id: String(row?.__index__ ?? index),
            } as IGVBaseFeature);
        }
    });

    return { features, sampled };
}

export { applyViewMargins, getLocationFieldsFromGenome, locationFromFieldValues };


