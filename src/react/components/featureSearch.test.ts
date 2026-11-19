import { describe, expect, it } from "vitest";
import { featureSearchMatcher, findMissingTerms, mergeTerms, parseFeatureList } from "./featureSearch";

describe("parseFeatureList", () => {
    it("splits on commas, semicolons, tabs and newlines", () => {
        expect(parseFeatureList("CD3E, CD4;CD8A\tMKI67\n\nPTPRC ")).toEqual(["CD3E", "CD4", "CD8A", "MKI67", "PTPRC"]);
    });
});

describe("mergeTerms", () => {
    it("ignores case-insensitive duplicates and keeps the first spelling", () => {
        expect(mergeTerms(["CD3E"], ["cd3e", "Cd4", "CD4"])).toEqual(["CD3E", "Cd4"]);
    });
});

describe("featureSearchMatcher", () => {
    it("matches everything when empty", () => {
        expect(featureSearchMatcher(undefined)).toBeNull();
        expect(featureSearchMatcher({ query: "  ", terms: [] })).toBeNull();
    });

    it("matches a substring of the query, ignoring case", () => {
        const m = featureSearchMatcher({ query: "cd" });
        expect(["CD3E", "PTPRC", "Cd4"].filter((v) => m?.(v))).toEqual(["CD3E", "Cd4"]);
    });

    it("matches terms exactly, ignoring case", () => {
        const m = featureSearchMatcher({ terms: ["cd3e", "MKI67"] });
        expect(["CD3E", "CD3EAP", "MKI67"].filter((v) => m?.(v))).toEqual(["CD3E", "MKI67"]);
    });

    it("narrows terms by the query", () => {
        const m = featureSearchMatcher({ terms: ["CD3E", "MKI67"], query: "ki" });
        expect(["CD3E", "MKI67"].filter((v) => m?.(v))).toEqual(["MKI67"]);
    });
});

describe("findMissingTerms", () => {
    it("lists terms not present in the column", () => {
        expect(findMissingTerms(["cd3e", "FOO", "MKI67"], ["CD3E", "MKI67"])).toEqual(["FOO"]);
    });
});
