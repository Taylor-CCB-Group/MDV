type JSONObject = { [key: string]: any };

/** Split a search query into lower-case words; all of them must match for something to be a hit. */
export function queryWords(query: string): string[] {
    return query.toLowerCase().split(/\s+/).filter(Boolean);
}

/** Whether every word of `query` appears (case-insensitively) somewhere in the given texts. An empty query matches. */
export function matchesQuery(query: string, ...texts: (string | undefined)[]): boolean {
    const words = queryWords(query);
    if (words.length === 0) return true;
    const haystack = texts.join(" ").toLowerCase();
    return words.every((w) => haystack.includes(w));
}

/**
 * Given a filter string and a json object, return a new json object that is a subset of the input.
 * A node is included if all (case insensitive) words of the filter appear somewhere in the node or its path;
 * a matching node is included with its whole subtree.
 */
export function filterJSON(obj: JSONObject, filter: string, ancestors: readonly string[] = []): JSONObject {
    const words = queryWords(filter);
    if (words.length === 0) return obj;

    function matchesFilter(value: string, path: string[]): boolean {
        const combinedString = path.concat(value).join(" ").toLowerCase();
        return words.every((w) => combinedString.includes(w));
    }

    function recursiveFilter(current: any, path: string[]): any {
        // return the whole subtree if the current node matches the filter
        if (matchesFilter("", path)) return current;
        if (typeof current === "string") {
            return matchesFilter(current, path) ? current : null;
        }
        if (typeof current === "object" && current !== null) {
            const filteredObj: JSONObject = {};
            let childMatch = false;
            for (const key in current) {
                const result = recursiveFilter(current[key], path.concat(key));
                if (result !== null) {
                    filteredObj[key] = result;
                    childMatch = true;
                }
            }
            if (childMatch) return filteredObj;
        }
        return null;
    }

    return recursiveFilter(obj, [...ancestors]) || {};
}

export type JsonSection = {
    /** Section title: the top-level key, or "General" for the collected primitive-valued keys. */
    title: string;
    json: unknown;
    /** For a group (see `splitJsonSections`): the sections nested inside it. `json` is unused for groups. */
    children?: JsonSection[];
};

/**
 * Split a debug JSON object into sections: each top-level key with an object/array value gets its own
 * section (so `buildInfo`, `datasources`, ... can be found by name), and primitive-valued keys are collected
 * into a "General" section.
 *
 * Keys listed in `groupKeys` become groups instead: a section whose children are the group value's own
 * entries (e.g. `project: { datasources, views, state }`).
 */
export function splitJsonSections(json: unknown, groupKeys: readonly string[] = []): JsonSection[] {
    if (typeof json !== "object" || json === null || Array.isArray(json)) {
        return [{ title: "JSON", json }];
    }
    const sections: JsonSection[] = [];
    const general: JSONObject = {};
    for (const [key, value] of Object.entries(json)) {
        if (typeof value !== "object" || value === null) general[key] = value;
        else if (groupKeys.includes(key) && !Array.isArray(value)) {
            sections.push({ title: key, json: value, children: splitJsonSections(value, []) });
        } else sections.push({ title: key, json: value });
    }
    if (Object.keys(general).length > 0) sections.unshift({ title: "General", json: general });
    return sections;
}

/**
 * Apply `query` to a JSON section. A matching title keeps the whole section; otherwise only matching
 * descendants (by key path or string value) are kept. Returns null if nothing in the section matches.
 */
export function filterJsonSection(
    section: JsonSection,
    query: string,
    /** Titles of the groups this section is nested in; they count as part of its path when matching. */
    ancestors: readonly string[] = [],
): JsonSection | null {
    if (queryWords(query).length === 0 || matchesQuery(query, ...ancestors, section.title)) return section;
    if (section.children) {
        const path = [...ancestors, section.title];
        const children = section.children.flatMap((child) => {
            const filtered = filterJsonSection(child, query, path);
            return filtered ? [filtered] : [];
        });
        return children.length > 0 ? { ...section, children } : null;
    }
    const filtered = filterJSON({ [section.title]: section.json }, query, ancestors)[section.title];
    return filtered === undefined ? null : { title: section.title, json: filtered };
}
