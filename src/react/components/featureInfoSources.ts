import z from "zod";
import { lookupMarker } from "./proteinAliases";

export const ORGANISMS = {
    human: { label: "Human", taxon: 9606 },
    mouse: { label: "Mouse", taxon: 10090 },
} as const;
export type Organism = keyof typeof ORGANISMS;

export type FeatureInfoLink = { label: string; href: string };

/** What a card shows, whichever source it came from. */
export type FeatureInfo = {
    title: string;
    /** Short identifiers line, e.g. gene symbol and accession. */
    subtitle?: string;
    description?: string;
    /** Extra one-line facts, e.g. subcellular location. */
    facts?: { label: string; value: string }[];
    /** Caveat about how the identifier was resolved. */
    note?: string;
    links: FeatureInfoLink[];
};

type SourceOptions = { organism: Organism; signal?: AbortSignal };

type SourceDefinition = {
    label: string;
    fetchInfo: (id: string, options: SourceOptions) => Promise<FeatureInfo>;
    /** Links offered when fetching fails. */
    fallbackLinks: (id: string, options: SourceOptions) => FeatureInfoLink[];
};

const genecardsLink = (id: string) => ({
    label: "GeneCards",
    href: `https://www.genecards.org/cgi-bin/carddisp.pl?gene=${encodeURIComponent(id)}`,
});

// Only the fields we show; the API returns much more.
const geneNetworkResponseSchema = z.object({
    gene: z.object({ id: z.string(), name: z.string(), description: z.string() }),
});

async function fetchGeneNetwork(id: string, { signal }: SourceOptions): Promise<FeatureInfo> {
    const response = await fetch(`https://www.genenetwork.nl/api/v1/gene/${encodeURIComponent(id)}`, { signal });
    if (!response.ok) {
        const errorData = await response.json().catch(() => null);
        throw new Error(errorData?.message || `HTTP error! status: ${response.status}`);
    }
    const { gene } = geneNetworkResponseSchema.parse(await response.json());
    return {
        title: gene.name,
        subtitle: gene.id,
        description: gene.description,
        links: [{ label: "GeneNetwork", href: `https://www.genenetwork.nl/gene/${gene.id}` }, genecardsLink(id)],
    };
}

const uniProtEntrySchema = z.object({
    primaryAccession: z.string(),
    proteinDescription: z
        .object({
            recommendedName: z.object({ fullName: z.object({ value: z.string() }) }).optional(),
            submissionNames: z.array(z.object({ fullName: z.object({ value: z.string() }) })).optional(),
        })
        .optional(),
    genes: z.array(z.object({ geneName: z.object({ value: z.string() }).optional() })).optional(),
    comments: z
        .array(
            z.object({
                commentType: z.string(),
                texts: z.array(z.object({ value: z.string() })).optional(),
                subcellularLocations: z
                    .array(z.object({ location: z.object({ value: z.string() }).optional() }))
                    .optional(),
            }),
        )
        .optional(),
});
const uniProtSearchSchema = z.object({ results: z.array(uniProtEntrySchema) });
type UniProtEntry = z.infer<typeof uniProtEntrySchema>;

const UNIPROT_FIELDS = "accession,gene_primary,protein_name,cc_function,cc_subcellular_location";

async function searchUniProt(query: string, signal?: AbortSignal) {
    const params = new URLSearchParams({ query, fields: UNIPROT_FIELDS, format: "json", size: "1" });
    const response = await fetch(`https://rest.uniprot.org/uniprotkb/search?${params}`, { signal });
    if (!response.ok) throw new Error(`UniProt error: status ${response.status}`);
    return uniProtSearchSchema.parse(await response.json()).results[0];
}

const uniProtSearchLink = (id: string, organism: Organism) => ({
    label: "Search UniProt",
    href: `https://www.uniprot.org/uniprotkb?query=${encodeURIComponent(`${id} AND organism_id:${ORGANISMS[organism].taxon}`)}`,
});

/** Strips characters that are UniProt query syntax, so an identifier can't change the query's meaning. */
const queryTerm = (s: string) => s.replace(/["():\\]/g, " ").trim();

function toFeatureInfo(entry: UniProtEntry, id: string, organism: Organism, note?: string): FeatureInfo {
    const gene = entry.genes?.[0]?.geneName?.value;
    const name =
        entry.proteinDescription?.recommendedName?.fullName.value ??
        entry.proteinDescription?.submissionNames?.[0]?.fullName.value ??
        gene ??
        entry.primaryAccession;
    const fn = entry.comments?.find((c) => c.commentType === "FUNCTION")?.texts?.[0]?.value;
    const locations = entry.comments
        ?.filter((c) => c.commentType === "SUBCELLULAR LOCATION")
        .flatMap((c) => c.subcellularLocations ?? [])
        .map((l) => l.location?.value)
        .filter((v): v is string => !!v);
    const links = [{ label: "UniProt", href: `https://www.uniprot.org/uniprotkb/${entry.primaryAccession}/entry` }];
    if (organism === "human" && gene) {
        links.push({ label: "Protein Atlas", href: `https://www.proteinatlas.org/search/${encodeURIComponent(gene)}` });
    }
    const idDiffers = gene && gene.toLowerCase() !== id.toLowerCase();
    return {
        title: name,
        subtitle: [idDiffers ? `${id} → ${gene}` : gene, entry.primaryAccession].filter(Boolean).join(" · "),
        description: fn,
        facts: locations?.length ? [{ label: "Location", value: [...new Set(locations)].join("; ") }] : undefined,
        note,
        links,
    };
}

async function fetchUniProt(id: string, { organism, signal }: SourceOptions): Promise<FeatureInfo> {
    const alias = lookupMarker(id);
    if (alias?.kind === "note") {
        return { title: id, description: alias.note, links: [uniProtSearchLink(id, organism)] };
    }
    const scope = `organism_id:${ORGANISMS[organism].taxon} AND reviewed:true`;
    const ensembl = /^ENS[A-Z]*G\d+/.exec(id)?.[0];
    if (ensembl) {
        // no organism scope: the Ensembl prefix already fixes the species (ENSG human, ENSMUSG mouse)
        const entry = await searchUniProt(`xref:ensembl-${ensembl} AND reviewed:true`, signal);
        if (entry) return toFeatureInfo(entry, id, organism);
    }
    const symbol = queryTerm(alias?.gene ?? id);
    if (symbol) {
        const entry = await searchUniProt(`gene_exact:"${symbol}" AND ${scope}`, signal);
        if (entry) return toFeatureInfo(entry, id, organism);
    }
    const text = queryTerm(id);
    const guess = text ? await searchUniProt(`(${text}) AND ${scope}`, signal) : undefined;
    if (guess) return toFeatureInfo(guess, id, organism, "Best guess from a UniProt text search");
    throw new Error(`No reviewed UniProt entry for "${id}"`);
}

export const FEATURE_INFO_SOURCES = {
    genenetwork: {
        label: "GeneNetwork (genes, human)",
        fetchInfo: fetchGeneNetwork,
        fallbackLinks: (id) => [genecardsLink(id)],
    },
    uniprot: {
        label: "UniProt (proteins)",
        fetchInfo: fetchUniProt,
        fallbackLinks: (id, { organism }) => [uniProtSearchLink(id, organism)],
    },
} satisfies Record<string, SourceDefinition>;
export type FeatureInfoSource = keyof typeof FEATURE_INFO_SOURCES;
