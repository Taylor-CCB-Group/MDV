import {
    type GenomeMetadata,
} from "../genomicLocationUtils";


export interface IGVBaseFeature {
    chr: string;
    start: number;
    end: number;
    id: string;
    color?: string;
    chr1?: string;
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


export function buildBaseFeatures(
    rows: Record<string, unknown>[],
    genomicInfo: GenomeMetadata,
    maxInitialFeatures: number,
): IGVBaseFeature[] {
  
    const maxRows = Math.min(rows.length, maxInitialFeatures);
    //nice and easy
    if (genomicInfo.type === "interval") {
        const map = genomicInfo.columns as Record<string, string>;
        const features: IGVBaseFeature[] =  new Array(maxRows);
        for (let i = 0; i < maxRows; i++) {
            const row = rows[i];
            features[i]= {
                chr: String(row[map.chr]),
                start: Number(row[map.start]),
                end: Number(row[map.end]),
                id: String((row.__index__ )),
            }
        }
        return features
    }
    else if (genomicInfo.type === "sv") {
        const map = genomicInfo.columns as  Record<string, string>;
        //don't know array size upfront because of potential for one or two features per row depending on svtype
        const features: IGVBaseFeature[] =  [];

        for (let i = 0; i < maxRows; i++) {
            const row = rows[i];
          
            const chr1 = String(row[map.chr1]);
            const pos1 = Number(row[map.pos1]);
            const pos2 = Number(row[map.pos2]);
            const chr2 = String(row[map.chr2]);
            if (!chr1 || !chr2 || !Number.isFinite(pos1) || !Number.isFinite(pos2)) {
                continue;
            }
            const id = String(row.__index__ );
            const svtype = String(row[map.svtype]);
            const length = Number(row[map.length]);
            const baseFeature = {
                id,
                chr2,
                pos1, 
                pos2,
                svtype,
                length,
             
            };
            if (svtype === "TRA" || svtype === "BND") {
                features.push({
                    ...baseFeature,
                    id,
                    chr: chr1,
                    start: pos1,
                    end: pos1,
                });
            
                if (typeof chr2 === "string") {
                    features.push({
                        ...baseFeature,
                        id,
                        chr: chr2,
                        start: pos2,
                        end: pos2,
                    });
                }
                continue;
            }
            if (svtype === "INS") {
                features.push({
                    ...baseFeature,
                    chr: chr1,
                    start: pos1,
                    end: pos2,
                });
                continue;
            }
            features.push({
                ...baseFeature,
                chr:chr1,
                start: Math.min(pos1, pos2),
                end: Math.max(pos1, pos2),
            });
            
        }
        return features;
       
    }
    return []; 
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



