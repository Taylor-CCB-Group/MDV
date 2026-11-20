import { parseDelimitedString } from "@/lib/utils";

/**
 * What the Gene Info search box filters the view by.
 * A row passes when it is one of `terms` (if any) and contains `query` (if any), ignoring case.
 */
export type FeatureSearch = {
    /** Free text typed into the box; a substring match. */
    query?: string;
    /** Identifiers pasted or entered one at a time; exact matches. */
    terms?: string[];
};

export function isSearchActive(search: FeatureSearch | undefined): search is FeatureSearch {
    return !!search && (!!search.query?.trim() || (search.terms?.length ?? 0) > 0);
}

/** Splits pasted text into identifiers (commas, semicolons or whitespace). */
export function parseFeatureList(text: string) {
    return parseDelimitedString(text.replaceAll(";", " "));
}

/** Adds `added` to `terms`, dropping case-insensitive duplicates and keeping the first spelling. */
export function mergeTerms(terms: readonly string[], added: readonly string[]) {
    const seen = new Set(terms.map((t) => t.toLowerCase()));
    const merged = [...terms];
    for (const t of added) {
        const key = t.toLowerCase();
        if (seen.has(key)) continue;
        seen.add(key);
        merged.push(t);
    }
    return merged;
}

/** Returns a predicate over identifier values, or `null` when the search matches everything. */
export function featureSearchMatcher(search: FeatureSearch | undefined): ((value: string) => boolean) | null {
    if (!isSearchActive(search)) return null;
    const query = search.query?.trim().toLowerCase() ?? "";
    const terms = search.terms?.length ? new Set(search.terms.map((t) => t.toLowerCase())) : null;
    return (value) => {
        const v = value.toLowerCase();
        if (terms && !terms.has(v)) return false;
        return query === "" || v.includes(query);
    };
}

/** The terms that match no value in the column, in the order given. */
export function findMissingTerms(terms: readonly string[], values: Iterable<string>) {
    if (terms.length === 0) return [];
    const present = new Set<string>();
    for (const v of values) present.add(v.toLowerCase());
    return terms.filter((t) => !present.has(t.toLowerCase()));
}
