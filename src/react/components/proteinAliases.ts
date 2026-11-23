/**
 * Antibody-panel marker names (CODEX, IMC, CosMx protein...) that UniProt can't resolve as gene
 * names, or resolves to the wrong protein by text search. Keys are normalised by {@link normaliseMarker}.
 */
const MARKER_GENES: Record<string, string> = {
    ki67: "MKI67",
    asma: "ACTA2",
    sma: "ACTA2",
    alphasma: "ACTA2",
    cd3: "CD3E",
    cd8: "CD8A",
    cd45: "PTPRC",
    cd45ra: "PTPRC",
    cd45ro: "PTPRC",
    cd20: "MS4A1",
    cd31: "PECAM1",
    cd11b: "ITGAM",
    cd11c: "ITGAX",
    cd56: "NCAM1",
    cd66b: "CEACAM8",
    hladr: "HLA-DRA",
    pdl1: "CD274",
    pd1: "PDCD1",
    ctla4: "CTLA4",
    ecadherin: "CDH1",
    ecad: "CDH1",
    ncadherin: "CDH2",
    epcam: "EPCAM",
    vimentin: "VIM",
    betacatenin: "CTNNB1",
    bcatenin: "CTNNB1",
    granzymeb: "GZMB",
    lyve1: "LYVE1",
    podoplanin: "PDPN",
    collageniv: "COL4A1",
    collagen4: "COL4A1",
    gfap: "GFAP",
    iba1: "AIF1",
    neun: "RBFOX3",
    mpo: "MPO",
    tryptase: "TPSAB1",
    synaptophysin: "SYP",
    chromogranina: "CHGA",
    er: "ESR1",
    pr: "PGR",
    her2: "ERBB2",
    p53: "TP53",
    ck7: "KRT7",
    ck8: "KRT8",
    ck18: "KRT18",
    ck19: "KRT19",
    ck20: "KRT20",
    ck5: "KRT5",
    ck14: "KRT14",
    cd16: "FCGR3A",
    cd107a: "LAMP1",
    cd366: "HAVCR2",
    tim3: "HAVCR2",
    b7h4: "VTCN1",
    pdl2: "PDCD1LG2",
    tbet: "TBX21",
    vegf: "VEGFA",
    collagentypei: "COL1A1",
    collagetypei: "COL1A1",
    collagen1: "COL1A1",
};

/** Markers that are not one protein, with what to say instead. */
const MARKER_NOTES: Record<string, string> = {
    cd57: "CD57 (HNK-1): a carbohydrate epitope carried by several molecules, not a single protein.",
    panck: "Pan-cytokeratin: an antibody cocktail against several keratins, not a single protein.",
    pancytokeratin: "Pan-cytokeratin: an antibody cocktail against several keratins, not a single protein.",
    pankeratin: "Pan-keratin: an antibody cocktail against several keratins, not a single protein.",
    hladrdqdp: "Pan HLA class II (HLA-DR, -DQ and -DP chains), not a single protein.",
    cytokeratin: "Cytokeratin: usually an antibody cocktail against several keratins, not a single protein.",
    dapi: "Nuclear stain (DNA), not a protein.",
    hoechst: "Nuclear stain (DNA), not a protein.",
    dna1: "DNA intercalator channel, not a protein.",
    dna2: "DNA intercalator channel, not a protein.",
};

/** Lower-cases and drops punctuation, spaces and a Greek alpha/beta prefix, so "α-SMA", "aSMA" and "SMA" meet. */
export function normaliseMarker(marker: string) {
    return marker
        .toLowerCase()
        .replace(/^α/, "alpha")
        .replace(/^β/, "beta")
        .replace(/[^a-z0-9]/g, "");
}

export type MarkerAlias = { kind: "gene"; gene: string } | { kind: "note"; note: string };

export function lookupMarker(marker: string): MarkerAlias | undefined {
    const key = normaliseMarker(marker);
    const note = MARKER_NOTES[key];
    if (note) return { kind: "note", note };
    const gene = MARKER_GENES[key];
    if (gene) return { kind: "gene", gene };
    return undefined;
}
